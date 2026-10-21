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

Run the frequency-constrained ridge and conservative component-reconstruction
example with:

```sh
julia --project examples/ridge_components.jl
```

It tracks two chirps in separately supplied search ranges, reconstructs them
with ±8 Hz bands, reports ridge, waveform, similarity, and
coefficient-conservation metrics, and writes `examples/ridge_components.png`.

```bash
julia --project=examples examples/damped_oscillators.jl
```

builds two chirped oscillators with different exponential damping times,
extracts and plots both ridges, reconstructs their analytic components, fits
their damping times, reports ridge/damping/reconstruction errors and
similarities, and writes `examples/damped_oscillators.png`.
An alternative output path can be supplied as the first argument.

For a noisy fade-out stress test in which one component disappears long
before the other, run:

```bash
julia --project=examples examples/damping_fadeout_stress.jl
```

The resulting `examples/damping_fadeout_stress.png` shows the detectability
boundary, ridge behavior after the short-lived mode vanishes, separated
reconstructions, reliability-aware damping fits, and damping bias as the fit
window is extended into noise.

To compare deterministic relative-amplitude termination with local WSST
prominence and persistent absorbing support, run:

```bash
julia --project=examples examples/mode_termination_two_components.jl
```

The example preserves the raw ridge and reconstruction for debugging, compares
2%, 5%, and 10% amplitude thresholds and their damping fits, tests prominence
thresholds independently, applies a combined operational cutoff, and writes
`examples/mode_termination_two_components.png`.

The nonexponential-decay and strong-neighbor leakage validations are run with:

```bash
julia --project=examples examples/mode_termination_validation.jl
```

This writes `examples/mode_termination_validation.png` and reports signed
stopping-time error and damping bias versus frequency separation.

Compare first-order, exact second-order, and prediction-window-pruned tracking
on identical synthetic time-frequency evidence with:

```sh
julia --project=examples examples/predictive_ridge_strategies.jl
```

The script reports ridge RMSE for every strategy and writes
`examples/predictive_ridge_strategies.png`. It deliberately includes a
stronger rapidly wandering competitor so smoothness is tested against known
truth rather than judged only by appearance.
