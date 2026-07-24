# Rule: Target-Centric CMake Architecture & Build Safety

*Activation Mode: Glob (`**/CMakeLists.txt`, `**/*.cmake`)*

## 1. Target-Centric Dependency Scope
- **No Legacy Globals:** Never use legacy global commands such as `include_directories()`, `link_directories()`, or `add_definitions()`.
- **Target Directives:** Scope all build properties to specific CMake targets using `target_include_directories()`, `target_link_libraries()`, `target_compile_definitions()`, and `target_compile_options()`.
- **Interface Visibility:** Explicitly tag target attributes with exact visibility scopes:
  - `PRIVATE`: Internal build dependencies not exposed in public header files.
  - `PUBLIC`: Dependencies required for both building the target and compiling consuming headers.
  - `INTERFACE`: Header-only library dependencies required by downstream consumers only.

```cmake
# GOOD: Target-scoped property definition
add_library(correlation_calculators STATIC ${CALCULATOR_SOURCES})
target_include_directories(correlation_calculators 
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(correlation_calculators
    PUBLIC
        correlation_core
    PRIVATE
        OpenMP::OpenMP_CXX
        TBB::tbb
)
```

## 2. Compile Commands & Language Standard
- **Compilation Database:** Always set `set(CMAKE_EXPORT_COMPILE_COMMANDS ON)` to produce `compile_commands.json` for `clang-tidy`, `clangd`, and static analysis tools.
- **Language Standard Discipline:** Require C++23 (`set(CMAKE_CXX_STANDARD 23)`, `set(CMAKE_CXX_STANDARD_REQUIRED ON)`, `set(CMAKE_CXX_EXTENSIONS OFF)`).
- **Out-of-Source Build Guard:** Prohibit in-source compilation (`PROJECT_SOURCE_DIR` == `PROJECT_BINARY_DIR`). Prevent pollutions of source trees.

## 3. Package & Dependency Fetching
- **Static Third-Party Packaging:** Default to static linking (`FetchContent` / static targets) for standalone tools (`nfd`, `Slint`, `voro++`) to prevent runtime dynamic library loading failures on host platforms.
- **Option Scoping:** Standardize feature toggles using `option(NAME "Description" DEFAULT_VALUE)` and propagate compile flags clean via target definitions.
