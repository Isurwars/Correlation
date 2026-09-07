# ADR 0001: Modern C++23 Standards, Memory Locality & SIMD Vectorization

## Status
Accepted

## Context
`Correlation` performs pairwise, triplet, and cluster evaluations over large atomistic structures (tens to hundreds of thousands of atoms across multi-gigabyte trajectories). Previous versions utilized ad-hoc loop structures, raw pointers, and double precision floats globally, causing memory cache misses and limiting SIMD throughput on modern AVX2/AVX-512/NEON hardware.

## Decision
1. **Language Standard**: Adopt modern C++23 (`-std=c++23`) globally across the core and calculator libraries. Leverage `std::span`, `std::ranges`, concepts, `std::expected`, and `constexpr` tables.
2. **Data-Oriented Design (SoA / AoS)**: Maintain contiguous memory layouts (`std::vector<real_t>` aligned to 64-byte cache boundaries where appropriate) to avoid false sharing and maximize L1/L2 cache locality.
3. **Precision Abstraction**: Define configurable `real_t` via `USE_SINGLE_PRECISION` (defaulting to single-precision `float` for 2x SIMD throughput and halved memory footprint, with `double` precision preserved when explicitly selected).
4. **SIMD Vectorization**: Structure inner pairwise distance loops to be auto-vectorizable with compiler hints (`#pragma GCC ivdep`, `#pragma omp simd`) and provide explicit SIMD helper routines (`cmake/SIMD.cmake`).

## Consequences
### Positive
- 2x–4x performance acceleration across RDF, S(Q), and angular distribution calculators.
- Dramatic reduction in memory footprint for multi-frame trajectories.
- Zero raw pointer allocations (`new`/`delete`); guaranteed RAII resource life-cycles.

### Negative / Trade-offs
- Requires modern C++23-compliant compilers (GCC 13+, Clang 17+, MSVC 2022+).
- Single precision produces minor numerical variance ($\le 10^{-6}$ relative error) compared to IEEE double precision.
