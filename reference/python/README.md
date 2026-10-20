# Python reference environment

This directory is exclusively for reproducible reference inputs and later
cross-language comparisons. `SynchroSqueezing.jl` never loads Python at runtime.

Create an isolated environment and generate deterministic input fixtures:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python generate_signals.py
python generate_reference.py
```

The reference package is pinned to `ssqueezepy==0.6.6`. Milestone 03A uses its
generalized-Morse wavelet with `gamma=3`, `beta=20`, matching the Julia
standard path. The staged fixture contains analytical input, raw CWT,
derivative, phase frequency and mask, SST coordinates and matrix, scalar
admissibility inverse, and two known-band component reconstructions. Generated
arrays are ignored by Git and can be reproduced from the pinned environment.

ContinuousWavelets.jl and ssqueezepy still choose scale bounds and periodic
versus reflected boundaries differently. Raw matrices are therefore retained
for inspection but are not asserted elementwise equal. The validation compares
the mathematical normalization and errors against the same analytical signals.
