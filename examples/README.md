# Examples

The numerical examples use only the package environment. The plotting example
has an isolated environment so CairoMakie remains an optional dependency. It
also uses the pinned `ssqueezepy==0.6.6` reference environment:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
reference/python/.venv/bin/python -m pip install -r reference/python/requirements.txt
julia --project=examples examples/wsst_three_oscillators.jl
```

The last command reads `three_oscillators.toml` from both Julia and Python.
The Python adapter converts the Julia inverse-scale bank to ssqueezepy's
dilation convention, and the example checks that the sample arrays and output
frequency grids agree. It writes `examples/wsst_three_oscillators.png`; an
alternative output path can be supplied as the first command-line argument.
Set `PYTHON` to select another Python environment containing the pinned
dependencies.
