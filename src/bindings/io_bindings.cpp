#include "io_bindings.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"
#include <filesystem>
#include <pybind11/pybind11.h>
#include <stdexcept>

namespace py = pybind11;
using namespace correlation::readers;
using namespace correlation::core;

void init_io(py::module_ &mod) {
  mod.def(
      "read",
      [](const std::string &filepath) -> Trajectory {
        const std::filesystem::path path(filepath);
        auto extension = path.extension().string();

        auto *reader = ReaderFactory::instance().getReaderForExtension({.extension = extension, .filename = filepath});
        if (!reader) {
          throw std::runtime_error("No reader found for file extension: " + extension);
        }

        return reader->readTrajectory(filepath, nullptr);
      },
      py::arg("filepath"),
      "Read an atomic trajectory or structure file into a Trajectory object.\n\n"
      "Automatically detects format from file extension (e.g. .car, .arc, .lammps, .xyz).\n\n"
      "Parameters\n----------\n"
      "filepath : str\n"
      "    Path to the input structure or trajectory file.\n\n"
      "Returns\n-------\n"
      "Trajectory\n"
      "    Parsed trajectory containing one or more simulation Cell frames.\n\n"
      "Raises\n------\n"
      "RuntimeError\n"
      "    If no reader is available for the given file extension or parsing fails.");
}
