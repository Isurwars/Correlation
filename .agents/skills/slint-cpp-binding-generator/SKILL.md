---
name: slint-cpp-binding-generator
description: Manages C++20 integrations with Slint UI, CMake target linkage, thread-safe data models, and event loop dispatching.
---

# Slint C++ Integration & Data Binding Standards

This skill directs high-performance C++20 integrations with Slint UI components, focusing on memory safety, CMake targets, and non-blocking asynchronous updates.

## Core Directives

1. **CMake Target Integration:** Build `.slint` markup into native C++20 targets via `slint_target_sources()` in CMake build scripts.
2. **Type Mapping:** Map Slint models to high-efficiency C++ containers using `slint::VectorModel<T>` and `slint::SharedString`.
3. **Event Loop Synchronization:** Always update UI state from background calculation threads (e.g. OpenMP analysis loops) using `slint::invoke_from_event_loop(...)`.
4. **RAII Handlers & Callbacks:** Bind Slint callbacks to C++ lambdas or member functions using weak references or RAII object lifetimes to avoid dangling memory traps.

## C++ Integration Pattern Example

```cpp
#include "app_window.h" // Generated header from Slint
#include <slint/slint.h>
#include <memory>
#include <thread>

void runApplication() {
    auto app = AppWindow::create();

    // Bind UI Callback to C++ Logic
    app->on_start_analysis([app_weak = std::weak_ptr<AppWindow>(app)](double cutoff) {
        std::thread([app_weak, cutoff]() {
            // Perform computation on worker thread
            double result_gr = 2.45; // ... calculation logic

            // Safe async UI update on main thread
            slint::invoke_from_event_loop([app_weak, result_gr]() {
                if (auto app = app_weak.lock()) {
                    app->set_calculation_status(slint::SharedString("Complete"));
                    app->set_result_amplitude(result_gr);
                }
            });
        }).detach();
    });

    app->run();
}
```

## CMake Setup Pattern

```cmake
find_package(Slint REQUIRED)

add_executable(correlation_gui src/main.cpp)
slint_target_sources(correlation_gui ui/app_window.slint)
target_link_libraries(correlation_gui PRIVATE Slint::Slint CorrelationCore)
target_compile_features(correlation_gui PRIVATE cxx_std_20)
```
