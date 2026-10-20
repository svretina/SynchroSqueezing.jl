"""Generate deterministic analytical signals for cross-language fixtures.

This script deliberately generates inputs only. Transform fixtures begin after
the Julia backend's conventions are selected and recorded.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import ssqueezepy


def sample_times(samples: int, sample_rate: float) -> np.ndarray:
    return np.arange(samples, dtype=np.float64) / sample_rate


def build_signals(samples: int, sample_rate: float) -> dict[str, np.ndarray]:
    time = sample_times(samples, sample_rate)
    tone = np.cos(2 * np.pi * 64.0 * time)
    chirp = np.cos(2 * np.pi * (32.0 * time + 0.5 * 40.0 * time**2))
    return {"time": time, "tone": tone, "linear_chirp": chirp}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("generated/signals.npz"))
    parser.add_argument("--samples", type=int, default=2048)
    parser.add_argument("--sample-rate", type=float, default=1024.0)
    args = parser.parse_args()

    if args.samples < 2:
        raise ValueError("samples must be at least two")
    if args.sample_rate <= 0:
        raise ValueError("sample-rate must be positive")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **build_signals(args.samples, args.sample_rate))
    metadata = {
        "sample_rate_hz": args.sample_rate,
        "samples": args.samples,
        "endpoint": "excluded",
        "numpy_version": np.__version__,
        "ssqueezepy_version": ssqueezepy.__version__,
    }
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
