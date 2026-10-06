# Rule: Target-Centric CMake Architecture & Build Safety

*Activation Mode: Glob (`**/CMakeLists.txt`, `**/*.cmake`, `CMakePresets.json`)*

## 1. Target-Centric Dependency Scope
- **No Legacy Globals:** Never introduce legacy global commands such as `link_directories()`. Minimize global flags; prefer target-scoped commands.
- **Target Directives:** Scope build properties to specific CMake targets using `target_include_directories()`, `target_link_libraries()`, `target_compile_definitions()`, and `target_compile_options()`.
- **Interface Visibility:** Explicitly tag target attributes with exact visibility scopes:
  - `PRIVATE`: Internal build dependencies not exposed in public header files.
  - `PUBLIC`: Dependencies required for both building the target and compiling consuming headers.
  - `INTERFACE`: Header-only library dependencies required by downstream consumers only.

```cmake
# Target-scoped property definition matching workspace architecture
target_include_directories(correlation_lib
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(correlation_lib
    PUBLIC
        TBB::tbb
    PRIVATE
        calculators_obj
        readers_obj
        writers_obj
)
```

## 2. Compile Commands & Language Standard
- **Compilation Database:** Always preserve `set(CMAKE_EXPORT_COMPILE_COMMANDS ON)` to produce `compile_commands.json` for `clang-tidy`, `clangd`, and static analysis tools.
- **Language Standard Discipline:** Require C++23 (`set(CMAKE_CXX_STANDARD 23)`, `set(CMAKE_CXX_STANDARD_REQUIRED ON)`, `set(CMAKE_CXX_EXTENSIONS OFF)`).
- **Out-of-Source Build Guard:** Prohibit in-source compilation (`PROJECT_SOURCE_DIR` == `PROJECT_BINARY_DIR`). Prevent pollutions of source trees.

## 3. Package & Dependency Fetching
- **Dependency Strategy:**
  - Standalone utility libraries (`nfd`, `voro++`, `miniz`, `CLI11`) are fetched via `FetchContent` when not provided by host.
  - Core system frameworks (`TBB`, `Slint`, `GTest`, `FFTW3`) are resolved via `find_package(... REQUIRED)`.
- **Option Scoping:** Standardize feature toggles using `option(NAME "Description" DEFAULT_VALUE)` and propagate compile flags cleanly via target definitions.
