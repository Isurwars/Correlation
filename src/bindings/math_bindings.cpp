/**
 * @file math_bindings.cpp
 * @brief Pybind11 bindings for correlation::math types.
 *
 * Exposes:
 *   - KernelType  (enum)
 */

#include "math_bindings.hpp"
#include "math/Smoothing.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace correlation::math;

void init_math(py::module_ &m) {
    // ------------------------------------------------------------------
    // KernelType enum
    // ------------------------------------------------------------------
    py::enum_<KernelType>(m, "KernelType",
        "Kernel type used for post-processing smoothing of histograms.")
        .value("Gaussian",  KernelType::Gaussian,
               "Gaussian (normal) kernel — smooth, infinite support.")
        .value("Bump",      KernelType::Bump,
               "Infinitely-smooth bump function with compact support.")
        .value("Triweight", KernelType::Triweight,
               "Triweight polynomial kernel with compact support.")
        .export_values();
}
