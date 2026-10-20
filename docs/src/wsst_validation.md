# WSST validation

The deterministic suite covers zero and constant inputs, real and analytic
complex tones, two separated tones, linear/quadratic/exponential chirps,
amplitude modulation, noise, odd/even lengths, Float32/Float64, small-amplitude
scaling, near-DC and near-Nyquist tones, and rejected two-sided complex input.

The tests separate four kinds of evidence:

1. Algebraic identities: accepted and deposited complex sums, linear first
   moments, and reference/fused agreement.
2. Physical correctness: tone/chirp frequencies in Hz at multiple sample
   rates.
3. Reconstruction: deposited-only, deposited plus independently accumulated
   residual, and DC-complete inverses.
4. Robustness limits: finite output near the represented band edges and in
   noise, without imposing interior-tone error targets there.

For a periodic 64 Hz Float64 tone, both nearest and linear full inverses are
required to have relative error below `1e-3`; current tests are substantially
tighter. A grid-refinement test uses a 37.3 Hz tone and verifies that reducing
output spacing reduces quantized ridge-frequency error. The backend response
is recomputed for every voices-per-octave configuration, so inversion does not
rely on an empirically chosen scalar.

## Python comparison

`reference/python/requirements.txt` pins NumPy and `ssqueezepy==0.6.6`.
`generate_reference.py` emits deterministic inputs and staged transform data
with versioned JSON metadata. Python remains an external reference environment;
it is never loaded by package code or Julia tests.

Milestone 03A uses matching generalized-Morse parameters and L¹ normalization
in both implementations. Scale bounds and padding still differ, so raw matrix
extents are not claimed to be elementwise interchangeable. Cross-language
checks validate the admissibility constant, analytical frequencies,
reconstruction, and known-component errors. The Julia randomized
fixed-coefficient test isolates phase and deposition from both CWT backends.

Run the Julia suite with `julia --project -e 'using Pkg; Pkg.test()'`, examples
with `julia --project examples/wsst_tone.jl`, and warmed performance
measurements from the benchmark environment described in
`benchmarks/README.md`.
