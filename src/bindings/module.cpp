#include "core_bindings.hpp"
#include "io_bindings.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(correlation, m) {
    m.doc() = "Correlation: Liquid and Amorphous Solid Analysis Tool Python Bindings";

    init_core(m);
    init_io(m);
}
