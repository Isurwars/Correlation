---
name: pybind11-interop-standards
description: Directs high-performance Python bindings creation, zero-copy buffer transfers, GIL management, and C++/Python memory safety.
---

# Pybind11 Interoperability & Python Packaging Standards

This skill governs the implementation, optimization, and maintenance of Python C++ extension bindings located in `src/bindings/`.

## 1. Zero-Copy NumPy Buffer Integration

1. **Direct Pointer Passing:** Prefer `py::array_t<T, py::array::c_style | py::array::forcecast>` to receive or export multi-dimensional array data without copying underlying memory.
2. **Buffer Protocol Interface:** Map C++ `Cell` and `Trajectory` coordinate buffers directly to Python using `py::buffer_info`.

```cpp
// GOOD: Zero-copy NumPy array exposure
m.def("compute_distances", [](py::array_t<double> coords) {
    py::buffer_info buf = coords.request();
    if (buf.ndim != 2 || buf.shape[1] != 3) {
        throw std::runtime_error("Input array must have shape (N, 3)");
    }
    
    std::span<const Vector3D> points(
        reinterpret_cast<const Vector3D*>(buf.ptr), 
        buf.shape[0]
    );
    return DistanceCalculator::compute(points);
});
```

---

## 2. Global Interpreter Lock (GIL) Management

1. **Release GIL During Heavy Computations:** Always instantiate `py::gil_scoped_release` before calling CPU-intensive OpenMP, TBB, or SYCL GPU calculations in C++.
2. **Re-acquire GIL on Callback:** If C++ code invokes a Python progress callback or logger, re-acquire the GIL using `py::gil_scoped_acquire`.

```cpp
// GOOD: Releasing GIL during compute phase
m.def("calculate_rdf", [](Trajectory& traj, const RDFNormalizationParams& params) {
    py::gil_scoped_release release;
    RDFCalculator calc;
    return calc.calculate(traj, params); // Multi-threaded TBB compute loop
});
```

---

## 3. Class Exposure & Memory Ownership

1. **Holders & Shared Pointers:** Use `std::shared_ptr<T>` as the holder type (`py::class_<T, std::shared_ptr<T>>`) when instances are co-owned by Python and C++ runtimes.
2. **Implicit Conversions & Type Safety:** Keep bindings strictly type-checked and raise informative Python exceptions (`std::invalid_argument` -> `ValueError`, `std::out_of_range` -> `IndexError`).
