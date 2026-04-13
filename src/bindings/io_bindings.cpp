#include "io_bindings.hpp"
#include "readers/ReaderFactory.hpp"
#include "core/Trajectory.hpp"
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <filesystem>

namespace py = pybind11;
using namespace correlation::readers;
using namespace correlation::core;

void init_io(py::module_ &m) {
    m.def("read", [](const std::string &filepath) -> Trajectory {
        std::filesystem::path p(filepath);
        auto extension = p.extension().string();
        
        auto* reader = ReaderFactory::instance().getReaderForExtension(extension);
        if (!reader) {
            throw std::runtime_error("No reader found for file extension: " + extension);
        }
        
        return reader->readTrajectory(filepath, nullptr);    }, "Read a trajectory from a file.");
}
