"""Generate staged ssqueezepy 0.6.6 reference data for Milestone 03."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import ssqueezepy
from ssqueezepy import Wavelet, cwt, icwt, issq_cwt, ssq_cwt
from ssqueezepy.utils import adm_ssq


def spectral_derivative(coefficients: np.ndarray, sample_rate: float) -> np.ndarray:
    samples = coefficients.shape[-1]
    angular = 2 * np.pi * np.fft.fftfreq(samples, d=1 / sample_rate)
    return np.fft.ifft(
        np.fft.fft(coefficients, axis=-1) * (1j * angular)[None, :],
        axis=-1,
    )


def synthetic_deposition_fixture() -> dict[str, np.ndarray]:
    centers = np.array([10.0, 20.0, 30.0])
    coefficients = np.array([[1 + 2j], [2 - 1j], [-0.5 + 0.25j]])
    frequencies = np.array([[10.0], [15.0], [30.0]])
    derivative = 1j * 2 * np.pi * frequencies * coefficients
    weights = np.array([0.2, 0.3, 0.5])
    expected = np.zeros((3, 1), dtype=np.complex128)
    expected[0, 0] += weights[0] * coefficients[0, 0]
    middle = weights[1] * coefficients[1, 0]
    expected[0, 0] += 0.5 * middle
    expected[1, 0] += 0.5 * middle
    expected[2, 0] += weights[2] * coefficients[2, 0]
    return {
        "synthetic_centers": centers,
        "synthetic_coefficients": coefficients,
        "synthetic_derivative": derivative,
        "synthetic_weights": weights,
        "synthetic_linear_expected": expected,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("generated/wsst_reference.npz"))
    parser.add_argument("--samples", type=int, default=1024)
    parser.add_argument("--sample-rate", type=float, default=512.0)
    args = parser.parse_args()
    if args.samples < 2 or args.sample_rate <= 0:
        raise ValueError("samples and sample-rate must be positive")

    time = np.arange(args.samples, dtype=np.float64) / args.sample_rate
    exact_frequency = 32.0 + 24.0 * time
    signal = np.cos(2 * np.pi * (32.0 * time + 12.0 * time**2))
    wavelet_config = ("gmw", {"gamma": 3, "beta": 20, "norm": "bandpass"})
    wavelet = Wavelet(wavelet_config)
    scale_policy = "log"
    coefficients, scales = cwt(
        signal, wavelet=wavelet, scales=scale_policy, nv=16, fs=args.sample_rate
    )
    derivative = spectral_derivative(coefficients, args.sample_rate)
    denominator = np.abs(coefficients) ** 2
    phase_frequency = np.full(denominator.shape, np.nan)
    valid = denominator > 1e-12 * denominator.max()
    phase_frequency[valid] = (
        np.imag(derivative[valid] * np.conj(coefficients[valid]))
        / (2 * np.pi * denominator[valid])
    )
    squeezed, _, ssq_frequencies, ssq_scales = ssq_cwt(
        signal, wavelet=wavelet, scales=scale_policy, nv=16, fs=args.sample_rate
    )
    reconstructed = icwt(coefficients, wavelet=wavelet, scales=scales)
    standard_reconstructed = issq_cwt(squeezed, wavelet=wavelet)

    tone1 = np.cos(2 * np.pi * 40.0 * time)
    tone2 = 0.7 * np.cos(2 * np.pi * 96.0 * time)
    two_tone_sst, _, two_tone_frequencies, _ = ssq_cwt(
        tone1 + tone2,
        wavelet=wavelet,
        scales=scale_policy,
        nv=16,
        fs=args.sample_rate,
    )
    centers = np.column_stack([
        np.full(args.samples, np.argmin(abs(two_tone_frequencies - frequency)))
        for frequency in (40.0, 96.0)
    ])
    half_widths = np.column_stack([
        np.full(args.samples, max(1, np.count_nonzero(abs(two_tone_frequencies - frequency) <= 4.0) // 2))
        for frequency in (40.0, 96.0)
    ])
    two_tone_components = issq_cwt(
        two_tone_sst,
        wavelet=wavelet,
        cc=centers,
        cw=half_widths,
    )
    interior = slice(args.samples // 8, -args.samples // 8)
    relative_error = lambda estimate, exact: float(
        np.linalg.norm(estimate - exact) / np.linalg.norm(exact)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        time=time,
        signal=signal,
        exact_frequency=exact_frequency,
        cwt=coefficients,
        scales=scales,
        derivative=derivative,
        phase_frequency=phase_frequency,
        phase_valid=valid,
        sst=squeezed,
        sst_frequencies=ssq_frequencies,
        sst_scales=ssq_scales,
        reconstructed=reconstructed,
        standard_reconstructed=standard_reconstructed,
        admissibility=np.array(adm_ssq(wavelet)),
        two_tone_component1=tone1,
        two_tone_component2=tone2,
        two_tone_sst=two_tone_sst,
        two_tone_frequencies=two_tone_frequencies,
        two_tone_reconstructed=two_tone_components,
        **synthetic_deposition_fixture(),
    )
    metadata = {
        "numpy_version": np.__version__,
        "ssqueezepy_version": ssqueezepy.__version__,
        "samples": args.samples,
        "sample_rate_hz": args.sample_rate,
        "wavelet": "gmw",
        "wavelet_parameters": {"gamma": 3, "beta": 20, "norm": "bandpass"},
        "scale_policy": scale_policy,
        "boundary_note": "ssqueezepy defaults; not elementwise aligned with periodic ContinuousWavelets.jl",
        "admissibility": float(adm_ssq(wavelet)),
        "standard_inverse_relative_l2_error": relative_error(
            standard_reconstructed, signal
        ),
        "two_tone_component_relative_l2_error": [
            relative_error(two_tone_components[0], tone1),
            relative_error(two_tone_components[1], tone2),
        ],
        "two_tone_component_interior_relative_l2_error": [
            relative_error(two_tone_components[0, interior], tone1[interior]),
            relative_error(two_tone_components[1, interior], tone2[interior]),
        ],
    }
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
