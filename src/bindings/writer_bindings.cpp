/**
 * @file writer_bindings.cpp
 * @brief Pybind11 bindings for correlation::writers types.
 *
 * Exposes:
 *   - BaseWriter      (abstract base)
 *   - CSVWriter
 *   - write_csv()     (free-function convenience wrapper)
 *   - get_writer()    (free function wrapping WriterFactory::instance())
 *   - list_writers()  (free function)
 *
 * Note: WriterFactory is a move-only singleton — it is not exposed as a
 * Python class. Instead, module-level free functions provide access to
 * its registered writers.
 */

#include "writer_bindings.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "writers/BaseWriter.hpp"
#include "writers/CSVWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>

namespace py = pybind11;
using namespace correlation::writers;
using namespace correlation::analysis;

void init_writers(py::module_ &m) {
    // ------------------------------------------------------------------
    // BaseWriter — abstract base (not directly instantiable from Python)
    // ------------------------------------------------------------------
    py::class_<BaseWriter>(m, "BaseWriter",
        "Abstract base class for all file-format writers.")
        .def("get_name",       &BaseWriter::getName,
             "Display name of this writer (e.g. 'CSV', 'HDF5').")
        .def("get_extensions", &BaseWriter::getExtensions,
             "List of supported file extensions (e.g. ['.csv']).")
        .def("write",
            [](const BaseWriter &w, const std::string &base_path,
               const DistributionFunctions &df, bool smoothing) {
                w.write(base_path, df, smoothing);
            },
            py::arg("base_path"), py::arg("df"), py::arg("smoothing") = false,
            "Write distribution function data to file(s).\n\n"
            "Parameters\n----------\n"
            "base_path : str\n"
            "    Base name for output files (without extension).\n"
            "df : DistributionFunctions\n"
            "    The analysis results to write.\n"
            "smoothing : bool, optional\n"
            "    If True, also write smoothed data. Default is False.");

    // ------------------------------------------------------------------
    // CSVWriter
    // ------------------------------------------------------------------
    py::class_<CSVWriter, BaseWriter>(m, "CSVWriter",
        "Writes distribution function histograms to CSV files.\n\n"
        "For each histogram (e.g. g(r)) a separate .csv file is created\n"
        "with all partials as columns.")
        .def(py::init<>())
        .def("write_all_csvs",
            [](const CSVWriter &w, const std::string &base_path,
               const DistributionFunctions &df, bool write_smoothed) {
                w.writeAllCSVs(base_path, df, write_smoothed);
            },
            py::arg("base_path"), py::arg("df"),
            py::arg("write_smoothed") = false,
            "Write all available histograms to individual CSV files.");

    // ------------------------------------------------------------------
    // WriterFactory access — exposed as module-level free functions
    // because WriterFactory holds unique_ptrs and cannot be copied.
    // ------------------------------------------------------------------
    m.def("get_writer",
        [](const std::string &name) -> BaseWriter * {
            auto *w = WriterFactory::instance().getWriter(name);
            if (!w)
                throw std::runtime_error("No writer registered with name: " + name);
            return w;
        },
        py::arg("name"),
        py::return_value_policy::reference,
        "Look up a file writer by name (e.g. 'CSV', 'HDF5').\n\n"
        "Returns a BaseWriter reference into the factory's singleton registry.\n\n"
        "Parameters\n----------\n"
        "name : str\n"
        "    Format name of the desired writer.\n\n"
        "Raises\n------\n"
        "RuntimeError\n"
        "    If no writer with that name is registered.");

    m.def("get_writer_for_extension",
        [](const std::string &ext) -> BaseWriter * {
            auto *w = WriterFactory::instance().getWriterForExtension(ext);
            if (!w)
                throw std::runtime_error("No writer found for extension: " + ext);
            return w;
        },
        py::arg("extension"),
        py::return_value_policy::reference,
        "Look up a writer by file extension (e.g. '.csv').");

    m.def("list_writers",
        []() {
            std::vector<std::string> names;
            for (const auto &w : WriterFactory::instance().getWriters())
                names.push_back(w->getName());
            return names;
        },
        "Return the names of all registered file writers.");

    // ------------------------------------------------------------------
    // Convenience free functions
    // ------------------------------------------------------------------
    m.def("write_csv",
        [](const std::string &base_path, const DistributionFunctions &df,
           bool write_smoothed) {
            CSVWriter w;
            w.writeAllCSVs(base_path, df, write_smoothed);
        },
        py::arg("base_path"), py::arg("df"),
        py::arg("write_smoothed") = false,
        "Convenience function: write all histograms in *df* to CSV files.\n\n"
        "Parameters\n----------\n"
        "base_path : str\n"
        "    Base name for output files (e.g. 'output/sample').\n"
        "df : DistributionFunctions\n"
        "    The analysis results to export.\n"
        "write_smoothed : bool, optional\n"
        "    If True, also write smoothed data files. Default is False.");
}
