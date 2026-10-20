# Python reference environment

This directory is exclusively for reproducible reference inputs and later
cross-language comparisons. `SynchroSqueezing.jl` never loads Python at runtime.

Create an isolated environment and generate deterministic input fixtures:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python generate_signals.py
```

The reference package is pinned to `ssqueezepy==0.6.6`. Generated files are
ignored by Git; canonical transform fixtures will only be committed after the
selected Julia backend’s wavelet, scale, padding, and normalization conventions
are matched.
