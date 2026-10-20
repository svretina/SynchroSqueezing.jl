"""Generate a matched ssqueezepy WSST for the three-oscillator example."""

from __future__ import annotations

import argparse
from pathlib import Path
import tomllib

import numpy as np
from ssqueezepy import Wavelet, ssq_cwt
from ssqueezepy.wavelets import center_frequency


def build_signal(config: dict) -> tuple[np.ndarray, np.ndarray]:
    samples = int(config["samples"])
    sample_rate = float(config["sample_rate_hz"])
    time = np.arange(samples, dtype=np.float64) / sample_rate

    oscillator1 = config["oscillator_1"]
    envelope = float(oscillator1["amplitude"]) * (
        1
        + float(oscillator1["modulation_depth"])
        * np.cos(2 * np.pi * float(oscillator1["modulation_frequency_hz"]) * time)
    )
    component1 = envelope * np.cos(
        2 * np.pi * float(oscillator1["frequency_hz"]) * time
    )

    components = [component1]
    for name in ("oscillator_2", "oscillator_3"):
        oscillator = config[name]
        amplitude = float(oscillator["amplitude"])
        initial = float(oscillator["initial_frequency_hz"])
        rate = float(oscillator["chirp_rate_hz_per_second"])
        phase_cycles = initial * time + 0.5 * rate * time**2
        components.append(amplitude * np.cos(2 * np.pi * phase_cycles))

    return time, np.sum(components, axis=0)


def relative_decibels(coefficients: np.ndarray) -> np.ndarray:
    relative = np.abs(coefficients) / np.max(np.abs(coefficients))
    return 20 * np.log10(np.maximum(relative, 1e-6))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.config.open("rb") as stream:
        config = tomllib.load(stream)

    _, signal = build_signal(config)
    wavelet_parameters = config["wavelet"]
    samples = int(config["samples"])
    voices = int(config["voices_per_octave"])
    wavelet = Wavelet(
        (
            "gmw",
            {
                "gamma": float(wavelet_parameters["gamma"]),
                "beta": float(wavelet_parameters["beta"]),
                "norm": "bandpass",
            },
        ),
        N=samples,
    )
    frequency_range = (
        float(config["frequency_min_hz"]),
        float(config["frequency_max_hz"]),
    )

    # ContinuousWavelets.jl uses inverse scales s_j = 2^(j / voices) and a
    # mother response peaking at its configured center frequency. Convert that
    # bank to ssqueezepy's dilation-scale convention so the sampled filters
    # cover the same nominal Fourier frequencies.
    maximum_scale_step = int(np.floor(voices * np.log2(max(2, samples // 2)))) - 1
    inverse_scales = 2 ** (np.arange(maximum_scale_step + 1) / voices)
    peak_angular_frequency = center_frequency(
        wavelet,
        kind="peak",
        N=samples,
    )
    center_frequency_parameter = float(wavelet_parameters["center_frequency"])
    scales = (
        peak_angular_frequency
        * samples
        / (2 * np.pi * center_frequency_parameter * inverse_scales)
    )[::-1]

    squeezed, _, frequencies, _ = ssq_cwt(
        signal,
        wavelet=wavelet,
        scales=scales,
        nv=voices,
        fs=float(config["sample_rate_hz"]),
        ssq_freqs="linear",
        padtype=None,
        squeezing="sum",
        maprange=frequency_range,
        gamma=0.0,
        astensor=False,
    )

    # ssqueezepy returns this configuration from high to low frequency.
    frequencies = np.asarray(frequencies[::-1], dtype=np.float64)
    squeezed = np.asarray(squeezed[::-1, :], dtype=np.complex128)
    expected_bins = int(config["frequency_bins"])
    if len(frequencies) != expected_bins:
        raise RuntimeError(
            f"expected {expected_bins} frequency bins, got {len(frequencies)}"
        )
    expected_frequencies = np.linspace(*frequency_range, expected_bins)
    if not np.allclose(frequencies, expected_frequencies, rtol=0, atol=1e-12):
        raise RuntimeError("ssqueezepy did not produce the requested frequency grid")

    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output / "signal.csv", signal, delimiter=",")
    np.savetxt(args.output / "frequencies.csv", frequencies, delimiter=",")
    np.savetxt(
        args.output / "wsst_decibels.csv",
        relative_decibels(squeezed),
        delimiter=",",
    )


if __name__ == "__main__":
    main()
