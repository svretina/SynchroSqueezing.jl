# Graph Report - .  (2026-10-20)

## Corpus Check
- 9 files · ~15,708 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 37 nodes · 39 edges · 8 communities detected
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]

## God Nodes (most connected - your core abstractions)
1. `BackendExperimentSupport` - 7 edges
2. `sample_times()` - 3 edges
3. `build_signals()` - 3 edges
4. `LinearAlgebra` - 3 edges
5. `compare_backends()` - 3 edges
6. `SynchroSqueezing` - 2 edges
7. `SynchroSqueezing` - 2 edges
8. `constant_tone()` - 2 edges
9. `linear_chirp()` - 2 edges
10. `sample_times()` - 2 edges

## Surprising Connections (you probably didn't know these)
- None detected - all connections are within the same source files.

## Communities

### Community 0 - "Community 0"
Cohesion: 0.29
Nodes (2): Statistics, BackendExperimentSupport

### Community 1 - "Community 1"
Cohesion: 0.47
Nodes (5): BenchmarkTools, ContinuousWavelet, compare_backends(), continuouswavelets_configuration(), trial_summary()

### Community 2 - "Community 2"
Cohesion: 0.4
Nodes (3): Documenter, SynchroSqueezing, Test

### Community 3 - "Community 3"
Cohesion: 0.6
Nodes (4): build_signals(), main(), Generate deterministic analytical signals for cross-language fixtures.  This scr, sample_times()

### Community 4 - "Community 4"
Cohesion: 0.4
Nodes (3): ContinuousWavelets, LinearAlgebra, Printf

### Community 5 - "Community 5"
Cohesion: 0.83
Nodes (3): constant_tone(), linear_chirp(), sample_times()

### Community 6 - "Community 6"
Cohesion: 0.67
Nodes (1): SynchroSqueezing

### Community 7 - "Community 7"
Cohesion: 1.0
Nodes (0):

## Knowledge Gaps
- **6 isolated node(s):** `Test`, `Documenter`, `Generate deterministic analytical signals for cross-language fixtures.  This scr`, `Statistics`, `BenchmarkTools` (+1 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 7`** (2 nodes): `main()`, `backend_comparison.jl`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `BackendExperimentSupport` connect `Community 0` to `Community 4`?**
  _High betweenness centrality (0.129) - this node is a cross-community bridge._
- **Why does `LinearAlgebra` connect `Community 4` to `Community 0`, `Community 1`?**
  _High betweenness centrality (0.117) - this node is a cross-community bridge._
- **Why does `ContinuousWavelets` connect `Community 4` to `Community 1`?**
  _High betweenness centrality (0.006) - this node is a cross-community bridge._
- **What connects `Test`, `Documenter`, `Generate deterministic analytical signals for cross-language fixtures.  This scr` to the rest of the system?**
  _6 weakly-connected nodes found - possible documentation gaps or missing edges._
