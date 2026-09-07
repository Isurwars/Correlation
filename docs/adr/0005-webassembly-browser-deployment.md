# ADR 0005: WebAssembly & Cross-Platform Browser Deployment via Emscripten

## Status
Accepted

## Context
Potential users and reviewers frequently require instant evaluation of scientific software without compiling native dependencies or installing local software suites. WebAssembly (WASM) enables executing modern C++ directly inside modern web browsers with near-native performance.

## Decision
1. **Toolchain**: Use Emscripten (`emcc`/`em++`) with CMake (`emcmake`) via the `BUILD_WASM=ON` CMake option.
2. **Embind Architecture**:
   - Bind C++ core classes (`Cell`, `Trajectory`, `Histogram`, `AnalysisSettings`, `DistributionFunctions`) and virtual file parsing (`readFromBuffer`) in `src/bindings/wasm_bindings.cpp`.
   - Expose typed views (`typed_memory_view`) to allow JavaScript/Canvas direct, zero-copy access to histogram arrays (`bins`, `partials`).
3. **Compilation Options**:
   - Single precision (`USE_SINGLE_PRECISION=ON`) to minimize WASM heap memory footprint.
   - SIMD128 enabled (`-msimd128`) for WebAssembly vector acceleration.
   - Modularized output (`-s MODULARIZE=1 -s EXPORT_NAME="createCorrelationModule"`).
   - Dynamic memory growth up to 2GB (`-s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=2GB`).
   - Virtual in-memory filesystem enabled (`-s FORCE_FILESYSTEM=1`) for seamless reader support.
4. **CI & Hosting**: Automate build and deployment to GitHub Pages via `.github/workflows/wasm.yml` leveraging `mymindstorm/setup-emsdk@v14`.

## Consequences
### Positive
- Zero-install, browser-based live demonstration and structural evaluation tool (`ui/wasm_app/`).
- Shared calculation algorithms between native desktop, Python bindings, and web runtimes.

### Negative / Trade-offs
- GUI-specific features (Slint desktop UI) and multi-threading (oneTBB pthreads) are excluded in the baseline WASM target to maintain broad browser compatibility without requiring COOP/COEP HTTP header constraints.
