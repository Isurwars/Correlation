#include "io_bindings.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"
#include "writers/structure/BaseStructureWriter.hpp"
#include "writers/structure/StructureConverter.hpp"
#include "writers/structure/StructureWriterFactory.hpp"

#include <filesystem>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;
using namespace correlation::readers;
using namespace correlation::core;
using namespace correlation::writers::structure;

namespace {

BaseStructureWriter *resolveWriter(const std::string &filepath, const std::string &format) {
  const auto &factory = StructureWriterFactory::instance();
  if (!format.empty()) {
    auto *writer = factory.getWriterForFormat(format);
    if (writer != nullptr) {
      return writer;
    }
    writer = factory.getWriterForExtension(format);
    if (writer != nullptr) {
      return writer;
    }
  }

  const std::filesystem::path path(filepath);
  std::string ext = path.extension().string();
  if (ext.empty()) {
    ext = path.filename().string();
  }
  return factory.getWriterForExtension(ext);
}

} // namespace

void init_io(py::module_ &mod) {
  mod.def(
      "read",
      [](const std::string &filepath) -> Trajectory {
        const std::filesystem::path path(filepath);
        const auto extension = path.extension().string();

        auto *reader = ReaderFactory::instance().getReaderForExtension(
            {.extension = extension, .filename = filepath});
        if (reader == nullptr) {
          throw std::runtime_error("No reader found for file extension: " + extension);
        }

        try {
          return reader->readTrajectory(filepath, nullptr);
        } catch (const std::exception &) {
          try {
            auto cell = reader->readStructure(filepath, nullptr);
            return Trajectory({std::move(cell)}, 1.0);
          } catch (...) {
            throw;
          }
        }
      },
      py::arg("filepath"),
      "Read an atomic trajectory or structure file into a Trajectory object.\n\n"
      "Automatically detects format from file extension (e.g. .car, .arc, .lammps, .xyz).\n\n"
      "Parameters\n----------\n"
      "filepath : str\n"
      "    Path to the input structure or trajectory file, or wildcard pattern.\n\n"
      "Returns\n-------\n"
      "Trajectory\n"
      "    Parsed trajectory containing one or more simulation Cell frames.\n\n"
      "Raises\n------\n"
      "RuntimeError\n"
      "    If no reader is available for the given file extension or parsing fails.");

  mod.def(
      "read",
      [](const std::vector<std::string> &filepaths) -> Trajectory {
        if (filepaths.empty()) {
          throw std::runtime_error("No files provided to read.");
        }
        const std::filesystem::path path(filepaths.front());
        const auto extension = path.extension().string();

        auto *reader = ReaderFactory::instance().getReaderForExtension(
            {.extension = extension, .filename = filepaths.front()});
        if (reader == nullptr) {
          throw std::runtime_error("No reader found for file extension: " + extension);
        }

        return reader->readTrajectory(filepaths, nullptr);
      },
      py::arg("filepaths"),
      "Read multiple atomic trajectory files into a single continuous Trajectory object.\n\n"
      "Parameters\n----------\n"
      "filepaths : list of str\n"
      "    List of paths to trajectory chunk files.\n\n"
      "Returns\n-------\n"
      "Trajectory\n"
      "    Parsed trajectory containing all frames concatenated in chronological order.");

  mod.def(
      "read_trajectory",
      [](const std::string &filepath) -> Trajectory {
        const std::filesystem::path path(filepath);
        const auto extension = path.extension().string();

        auto *reader = ReaderFactory::instance().getReaderForExtension(
            {.extension = extension, .filename = filepath});
        if (reader == nullptr) {
          throw std::runtime_error("No reader found for file extension: " + extension);
        }

        return reader->readTrajectory(filepath, nullptr);
      },
      py::arg("filepath"), "Read a trajectory file or wildcard pattern into a Trajectory object.");

  mod.def(
      "read_trajectory",
      [](const std::vector<std::string> &filepaths) -> Trajectory {
        if (filepaths.empty()) {
          throw std::runtime_error("No files provided to read.");
        }
        const std::filesystem::path path(filepaths.front());
        const auto extension = path.extension().string();

        auto *reader = ReaderFactory::instance().getReaderForExtension(
            {.extension = extension, .filename = filepaths.front()});
        if (reader == nullptr) {
          throw std::runtime_error("No reader found for file extension: " + extension);
        }

        return reader->readTrajectory(filepaths, nullptr);
      },
      py::arg("filepaths"),
      "Read multiple trajectory files into a single continuous Trajectory object.");

  mod.def(
      "read_structure",
      [](const std::string &filepath) -> Cell {
        const std::filesystem::path path(filepath);
        const auto extension = path.extension().string();

        auto *reader = ReaderFactory::instance().getReaderForExtension(
            {.extension = extension, .filename = filepath});
        if (reader == nullptr) {
          throw std::runtime_error("No reader found for file extension: " + extension);
        }

        return reader->readStructure(filepath, nullptr);
      },
      py::arg("filepath"),
      "Read an atomic structure file directly into a Cell object.\n\n"
      "Parameters\n----------\n"
      "filepath : str\n"
      "    Path to the input structure file (e.g. .car, POSCAR, .cif, .cell).\n\n"
      "Returns\n-------\n"
      "Cell\n"
      "    Parsed simulation cell.");

  mod.def(
      "write",
      [](const Cell &cell, const std::string &filepath, const std::string &format, bool fractional,
         const std::string &title) {
        const auto *const writer = resolveWriter(filepath, format);
        if (writer == nullptr) {
          throw std::runtime_error("No structure writer found for filepath: " + filepath);
        }
        writer->write(filepath, cell, {.fractional = fractional, .title = title});
      },
      py::arg("cell"), py::arg("filepath"), py::arg("format") = "", py::arg("fractional") = true,
      py::arg("title") = "",
      "Write a Cell object to an atomic structure file.\n\n"
      "Parameters\n----------\n"
      "cell : Cell\n"
      "    The simulation cell to serialize.\n"
      "filepath : str\n"
      "    Destination file path.\n"
      "format : str, optional\n"
      "    Target format identifier (e.g. 'poscar', 'cell', 'cif', 'car').\n"
      "fractional : bool, default True\n"
      "    Whether to use fractional (direct) coordinates if supported.\n"
      "title : str, optional\n"
      "    Optional title or header comment.");

  mod.def(
      "convert",
      [](const std::string &input_file, const std::string &output_file,
         const std::string &from_format, const std::string &to_format, bool fractional) {
        correlation::writers::structure::StructureConvertOptions const opts{
            .from_format = from_format,
            .to_format = to_format,
            .fractional = fractional,
        };
        correlation::writers::structure::convertStructure(input_file, output_file, opts);
      },
      py::arg("input_file"), py::arg("output_file"), py::arg("from_format") = "",
      py::arg("to_format") = "", py::arg("fractional") = true,
      "Convert directly between atomic structure and cell formats.\n\n"
      "Parameters\n----------\n"
      "input_file : str\n"
      "    Source structure file path.\n"
      "output_file : str\n"
      "    Destination structure file path.\n"
      "from_format : str, optional\n"
      "    Source format override.\n"
      "to_format : str, optional\n"
      "    Destination format override.\n"
      "fractional : bool, default True\n"
      "    Whether to format coordinates in fractional space.");
}
