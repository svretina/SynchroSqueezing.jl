# SynchroSqueezing.jl roadmap

This roadmap turns `instructions.md` into sequential, independently verifiable milestones.

The file contains two proposals. The latter, **“Revised architecture: Synchrosqueezing.jl using the Julia ecosystem,”** controls this roadmap whenever it conflicts with the earlier proposal. In particular, production code should reuse mature Julia CWT, STFT, FFT, and signal-processing implementations where they meet the scientific requirements. Synchrosqueezing, derivative integration, ridge tracking, and SST reconstruction remain this package’s responsibility. The repository and Julia package are named `SynchroSqueezing.jl` / `SynchroSqueezing`.

| Milestone | Outcome | Gate |
| --- | --- | --- |
| [00](00_foundations_and_dependency_audit.md) | Reproducible package foundation and dependency decision | Candidate backends are documented and the package/test/docs skeleton works. |
| [01](01_cwt_backend_prototypes.md) | Measured CWT backend choice | CWT behavior, inversion, layouts, and repeated-workload trade-offs are known. |
| [02](02_derivative_and_backend_interface.md) | SST-ready backend adapter | The physical-time CWT derivative and metadata contract are validated. |
| [03](03_first_order_wsst.md) | Minimum viable first-order WSST | Phase estimation, complex reassignment, and inverse WSST pass analytical tests. |
| [04](04_ridges_and_component_recovery.md) | Scientific mode recovery | Ridges and components are recoverable with explicit residual diagnostics. |
| [05](05_cpu_optimization_and_benchmarks.md) | Reproducible performant CPU pipeline | Profiled, low-allocation, race-free performance results are published. |
| [06](06_fsst_sst2_multiridge_and_bns.md) | Advanced methods and BNS validation | FSST/SST2 and crossing-aware multiridge behavior are analytically characterized. |

## Global rules

- Julia is the complete production runtime. Python is limited to reference fixtures and comparisons with a pinned `ssqueezepy` version.
- Establish and document Fourier signs, sampling units, normalization, scale/frequency grids, coefficient orientation, reassignment weights, and inverse constants before claiming numerical parity.
- Analytical tests are required in addition to cross-language comparisons.
- Never silently change mathematical conventions or tune normalizations merely to match an external output.
- Implement a scalar reference before optimizing; report measured performance fairly, including results where Python is faster.
- Keep optional packages optional unless the core runtime demonstrably needs them.
