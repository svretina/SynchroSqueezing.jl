#!/usr/bin/env python3
"""Pinned ssqueezepy half of the strict cross-language benchmark.

The Julia driver writes raw, typed fixtures.  File I/O and process startup are
outside every timed region; this program only times native Python calls.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import time
import tracemalloc
from pathlib import Path

import numba
import numpy as np
import scipy
import ssqueezepy
from ssqueezepy import Wavelet, cwt, issq_cwt
from ssqueezepy._ssq_cwt import phase_cwt
from ssqueezepy.ssqueezing import ssqueeze
from ssqueezepy.utils import adm_ssq


def load_array(root: Path, spec: dict, name: str) -> np.ndarray:
    item = spec["arrays"][name]
    dtype = np.dtype(item["dtype"])
    array = np.fromfile(root / item["file"], dtype=dtype)
    return array.reshape(item["shape"], order="F")


def error_metrics(actual: np.ndarray, reference: np.ndarray) -> dict:
    delta = np.asarray(actual) - np.asarray(reference)
    denominator = np.linalg.norm(reference.ravel())
    return {
        "relative_l2": float(np.linalg.norm(delta.ravel()) / max(denominator, np.finfo(float).tiny)),
        "maximum_absolute": float(np.max(np.abs(delta))),
    }


def timed(callable_, samples: int) -> tuple[dict, object]:
    first_start = time.perf_counter_ns()
    value = callable_()
    first_ms = (time.perf_counter_ns() - first_start) / 1e6
    callable_()  # warm Numba, FFT, and internal caches
    values = []
    tracemalloc.start()
    for _ in range(samples):
        start = time.perf_counter_ns()
        value = callable_()
        values.append((time.perf_counter_ns() - start) / 1e6)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    values.sort()
    return ({
        "first_call_ms": first_ms,
        "minimum_ms": min(values),
        "median_ms": statistics.median(values),
        "q25_ms": float(np.quantile(values, 0.25)),
        "q75_ms": float(np.quantile(values, 0.75)),
        "samples": len(values),
        "python_tracemalloc_peak_bytes": peak,
    }, value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=12)
    args = parser.parse_args()
    root = args.fixture
    spec = json.loads((root / "manifest.json").read_text())

    x = load_array(root, spec, "signal")
    julia_scales = load_array(root, spec, "scales")
    frequencies = load_array(root, spec, "frequencies")
    reference_wx = load_array(root, spec, "cwt")
    reference_dwx = load_array(root, spec, "derivative")
    reference_phase = load_array(root, spec, "phase")
    reference_tx = load_array(root, spec, "wsst")
    reference_inverse = load_array(root, spec, "inverse")

    n = int(spec["configuration"]["samples"])
    fs = float(spec["configuration"]["sample_rate_hz"])
    q = int(spec["configuration"]["voices_per_octave"])
    dtype_name = spec["configuration"]["precision"]
    real_dtype = np.float32 if dtype_name == "Float32" else np.float64
    complex_dtype = np.complex64 if real_dtype is np.float32 else np.complex128
    x = x.astype(real_dtype, copy=False)

    # ContinuousWavelets.jl's periodic coordinate ends at N/2 + 1 while
    # containing N/2 + 1 samples.  This endpoint factor is part of the sampled
    # filter bank and is therefore included in the exact scale adapter.
    gamma_parameter, beta_parameter = 3.0, 20.0
    peak = (beta_parameter / gamma_parameter) ** (1 / gamma_parameter)
    endpoint_factor = (n // 2 + 1) / (n // 2)
    mapped_descending = peak * n * endpoint_factor / (2 * np.pi * julia_scales)
    python_scales = mapped_descending[::-1].astype(real_dtype)
    wavelet = Wavelet(
        ("gmw", {"gamma": gamma_parameter, "beta": beta_parameter, "norm": "bandpass"}),
        N=n,
        dtype="float32" if real_dtype is np.float32 else "float64",
    )

    def cwt_only():
        wx, _ = cwt(x, wavelet=wavelet, scales=python_scales, fs=fs,
                    padtype=None, derivative=False, vectorized=True,
                    astensor=False, cache_wavelet=True, nan_checks=False)
        return np.asarray(wx)[::-1].astype(complex_dtype, copy=False)

    def cwt_derivative():
        wx, scales, dwx = cwt(x, wavelet=wavelet, scales=python_scales, fs=fs,
                              padtype=None, derivative=True, vectorized=True,
                              astensor=False, cache_wavelet=True, nan_checks=False)
        return (np.asarray(wx)[::-1].astype(complex_dtype, copy=False),
                np.asarray(dwx)[::-1].astype(complex_dtype, copy=False), scales)

    cwt_timing, wx = timed(cwt_only, args.samples)
    pair_timing, pair = timed(cwt_derivative, args.samples)
    wx_pair, dwx, returned_scales = pair
    threshold = float(spec["configuration"]["absolute_threshold"])

    def phase_only():
        return np.asarray(phase_cwt(wx_pair, dwx, "trig", threshold))

    phase_timing, phase = timed(phase_only, args.samples)
    valid = np.isfinite(reference_phase) & np.isfinite(phase)
    phase_error = error_metrics(phase[valid], reference_phase[valid])
    phase_error["compared_coefficients"] = int(valid.sum())

    # ssqueezepy requires its native ascending scale order.  Its reassignment
    # uses a constant log-scale weight, whereas Julia deliberately uses
    # trapezoidal half-weights at the two finite-bank endpoints.  We execute and
    # audit it, but the equivalence gate decides whether timings are comparable.
    native_wx = wx_pair[::-1]
    native_dwx = dwx[::-1]
    maprange = (float(frequencies[0]), float(frequencies[-1]))

    def reassignment_only():
        tx, returned_frequencies = ssqueeze(
            native_wx.copy(), w=None, ssq_freqs=frequencies.copy(),
            scales=python_scales, fs=fs, gamma=threshold, dWx=native_dwx,
            transform="cwt", squeezing="sum", maprange=maprange,
            was_padded=False, flipud=False,
        )
        # With flipud=false the kernel rows follow the supplied ascending
        # frequency array, although ssqueezepy reverses only the returned
        # coordinate vector. Keep the kernel rows in their actual order.
        return np.asarray(tx), np.asarray(returned_frequencies)

    reassignment_timing, reassigned = timed(reassignment_only, args.samples)
    tx, returned_frequencies = reassigned

    def complete_wsst():
        wx_native, _, dwx_native = cwt(
            x, wavelet=wavelet, scales=python_scales, fs=fs, padtype=None,
            derivative=True, vectorized=True, astensor=False,
            cache_wavelet=True, nan_checks=False,
        )
        tx_native, _ = ssqueeze(
            wx_native, w=None, ssq_freqs=frequencies.copy(),
            scales=python_scales, fs=fs, gamma=threshold, dWx=dwx_native,
            transform="cwt", squeezing="sum", maprange=maprange,
            was_padded=False, flipud=False,
        )
        return np.asarray(tx_native)

    complete_timing, complete_tx = timed(complete_wsst, max(4, args.samples // 2))
    admissibility = float(spec["configuration"]["admissibility"])
    python_admissibility = float(adm_ssq(wavelet))

    def inverse_only():
        return issq_cwt(reference_tx.copy(), wavelet=wavelet)

    inverse_timing, inverse = timed(inverse_only, args.samples)
    mask = load_array(root, spec, "mask").astype(bool)

    def band_inverse():
        return (2 / python_admissibility) * np.real(np.where(mask, reference_tx, 0).sum(axis=0))

    band_timing, band = timed(band_inverse, args.samples)
    reference_band = load_array(root, spec, "band_inverse")

    tolerances = spec["tolerances"]
    checks = {
        "cwt": error_metrics(wx, reference_wx),
        "cwt_from_combined_call": error_metrics(wx_pair, reference_wx),
        "derivative": error_metrics(dwx, reference_dwx),
        "phase": phase_error,
        "reassignment": error_metrics(tx, reference_tx),
        "complete_wsst": error_metrics(complete_tx, reference_tx),
        "standard_inverse": error_metrics(inverse, reference_inverse),
        "fixed_band_inverse": error_metrics(band, reference_band),
    }
    checks["reassignment"]["orientation_audit"] = {
        "current": checks["reassignment"]["relative_l2"],
        "row_reversed": error_metrics(tx[::-1], reference_tx)["relative_l2"],
        "complex_conjugated": error_metrics(np.conj(tx), reference_tx)["relative_l2"],
        "reference_peak_rows": [int(i) for i in np.argsort(np.linalg.norm(reference_tx, axis=1))[-3:]],
        "python_peak_rows": [int(i) for i in np.argsort(np.linalg.norm(tx, axis=1))[-3:]],
        "returned_frequencies_ascending": bool(np.all(np.diff(returned_frequencies) > 0)),
    }
    tolerance_names = {"cwt_from_combined_call": "cwt"}
    for name, check in checks.items():
        tolerance = tolerances[tolerance_names.get(name, name)]
        check["passed"] = (check["relative_l2"] <= tolerance["relative_l2"] and
                           check["maximum_absolute"] <= tolerance["maximum_absolute"])

    result = {
        "schema_version": 1,
        "implementation": "Python ssqueezepy",
        "versions": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "numba": numba.__version__,
            "ssqueezepy": ssqueezepy.__version__,
        },
        "environment": {
            "platform": platform.platform(), "processor": platform.processor(),
            "pid": os.getpid(), "numba_threads": numba.get_num_threads(),
            "thread_environment": {key: os.environ.get(key) for key in (
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMBA_NUM_THREADS", "SSQ_GPU")},
        },
        "configuration": spec["configuration"],
        "scale_adapter": {
            "formula": "a_py = omega_peak * N * ((N/2+1)/(N/2)) / (2*pi*s_julia)",
            "endpoint_factor": endpoint_factor,
            "python_scale_order": "ascending dilation; outputs reversed for ascending frequency",
            "nyquist_note": "ssqueezepy halves its analytic Nyquist filter bin; the deterministic fixture has zero Nyquist input.",
        },
        "admissibility": {
            "julia": admissibility,
            "python": python_admissibility,
            "absolute_difference": abs(admissibility - python_admissibility),
        },
        "equivalence": checks,
        "timings": {
            "cwt": cwt_timing,
            "cwt_and_derivative": pair_timing,
            "phase_estimation": phase_timing,
            "reassignment": reassignment_timing,
            "complete_first_order_wsst": complete_timing,
            "standard_inverse_wsst": inverse_timing,
            "fixed_band_reconstruction": band_timing,
        },
        "strict_comparability": {
            "cwt": checks["cwt"]["passed"],
            "cwt_and_derivative": checks["derivative"]["passed"],
            "phase_estimation": checks["phase"]["passed"],
            "reassignment": checks["reassignment"]["passed"],
            "complete_first_order_wsst": checks["complete_wsst"]["passed"],
            "standard_inverse_wsst": checks["standard_inverse"]["passed"],
            "fixed_band_reconstruction": checks["fixed_band_inverse"]["passed"],
            "ridge_tracking": False,
        },
        "unsupported": {
            "analytic_complex_inverse": "ssqueezepy issq_cwt exposes real synthesis only",
            "ridge_tracking": "objectives and transition constraints differ",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
