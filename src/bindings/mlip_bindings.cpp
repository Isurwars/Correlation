/**
 * @file mlip_bindings.cpp
 * @brief Pybind11 bindings for machine learning interatomic potentials (MLIP) and GNN graph builder.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip_bindings.hpp"
#include "calculators/MLIPCalculator.hpp"
#include "core/Cell.hpp"
#include "mlip/MLIPInterface.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using correlation::real_t;
using namespace correlation::mlip;
using namespace correlation::calculators;

void init_mlip(py::module_ &mod) {
  // ------------------------------------------------------------------
  // PeriodicGraphData
  // ------------------------------------------------------------------
  py::class_<PeriodicGraphData>(mod, "PeriodicGraphData",
                                "Container for periodic neighbor graph tensors ready for GNN model inference.")
      .def(py::init<>())
      .def_readonly("atom_count", &PeriodicGraphData::atom_count, "Total number of atoms N.")
      .def_readonly("edge_count", &PeriodicGraphData::edge_count, "Total number of directed edges E.")
      .def_property_readonly(
          "positions",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            if (graph.positions_flat.empty() || graph.atom_count == 0) {
              return py::array_t<real_t>();
            }
            py::ssize_t rows = static_cast<py::ssize_t>(graph.atom_count);
            py::ssize_t cols = 3;
            py::ssize_t stride_row = static_cast<py::ssize_t>(3 * sizeof(real_t));
            py::ssize_t stride_col = static_cast<py::ssize_t>(sizeof(real_t));
            return py::array_t<real_t>({rows, cols}, {stride_row, stride_col}, graph.positions_flat.data(), obj);
          },
          "Zero-copy access to atomic Cartesian coordinates as a (N, 3) NumPy array.")
      .def_property_readonly(
          "atomic_numbers",
          [](py::object &obj) -> py::array_t<int64_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            if (graph.atomic_numbers.empty() || graph.atom_count == 0) {
              return py::array_t<int64_t>();
            }
            py::ssize_t count = static_cast<py::ssize_t>(graph.atom_count);
            py::ssize_t stride = static_cast<py::ssize_t>(sizeof(int64_t));
            return py::array_t<int64_t>({count}, {stride}, graph.atomic_numbers.data(), obj);
          },
          "Zero-copy access to atomic numbers (Z) as a (N,) NumPy array.")
      .def_property_readonly(
          "edge_index",
          [](py::object &obj) -> py::array_t<int64_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            if (graph.edge_index_flat.empty() || graph.edge_count == 0) {
              return py::array_t<int64_t>();
            }
            py::ssize_t rows = 2;
            py::ssize_t cols = static_cast<py::ssize_t>(graph.edge_count);
            py::ssize_t stride_row = static_cast<py::ssize_t>(graph.edge_count * sizeof(int64_t));
            py::ssize_t stride_col = static_cast<py::ssize_t>(sizeof(int64_t));
            return py::array_t<int64_t>({rows, cols}, {stride_row, stride_col}, graph.edge_index_flat.data(), obj);
          },
          "Zero-copy access to directed edge indices (COO format) as a (2, E) NumPy array.")
      .def_property_readonly(
          "edge_shifts",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            if (graph.edge_shifts_flat.empty() || graph.edge_count == 0) {
              return py::array_t<real_t>();
            }
            py::ssize_t rows = static_cast<py::ssize_t>(graph.edge_count);
            py::ssize_t cols = 3;
            py::ssize_t stride_row = static_cast<py::ssize_t>(3 * sizeof(real_t));
            py::ssize_t stride_col = static_cast<py::ssize_t>(sizeof(real_t));
            return py::array_t<real_t>({rows, cols}, {stride_row, stride_col}, graph.edge_shifts_flat.data(), obj);
          },
          "Zero-copy access to periodic cell integer shift vectors as a (E, 3) NumPy array.")
      .def_property_readonly(
          "edge_vectors",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            if (graph.edge_vectors_flat.empty() || graph.edge_count == 0) {
              return py::array_t<real_t>();
            }
            py::ssize_t rows = static_cast<py::ssize_t>(graph.edge_count);
            py::ssize_t cols = 3;
            py::ssize_t stride_row = static_cast<py::ssize_t>(3 * sizeof(real_t));
            py::ssize_t stride_col = static_cast<py::ssize_t>(sizeof(real_t));
            return py::array_t<real_t>({rows, cols}, {stride_row, stride_col}, graph.edge_vectors_flat.data(), obj);
          },
          "Zero-copy access to Cartesian displacement vectors r_ij as a (E, 3) NumPy array.")
      .def_property_readonly(
          "edge_distances",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            if (graph.edge_distances.empty() || graph.edge_count == 0) {
              return py::array_t<real_t>();
            }
            py::ssize_t count = static_cast<py::ssize_t>(graph.edge_count);
            py::ssize_t stride = static_cast<py::ssize_t>(sizeof(real_t));
            return py::array_t<real_t>({count}, {stride}, graph.edge_distances.data(), obj);
          },
          "Zero-copy access to Euclidean edge distances ||r_ij|| as a (E,) NumPy array.")
      .def_property_readonly(
          "cell",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph = obj.cast<const PeriodicGraphData &>();
            py::ssize_t rows = 3;
            py::ssize_t cols = 3;
            py::ssize_t stride_row = static_cast<py::ssize_t>(3 * sizeof(real_t));
            py::ssize_t stride_col = static_cast<py::ssize_t>(sizeof(real_t));
            return py::array_t<real_t>({rows, cols}, {stride_row, stride_col}, graph.cell_flat.data(), obj);
          },
          "Zero-copy access to lattice matrix as a (3, 3) NumPy array.");

  // ------------------------------------------------------------------
  // GaussianRBFConfig
  // ------------------------------------------------------------------
  py::class_<GaussianRBFConfig>(mod, "GaussianRBFConfig",
                                "Configuration parameters for Gaussian radial basis function expansion.")
      .def(py::init<real_t, real_t, size_t>(), py::arg("start") = static_cast<real_t>(0.0),
           py::arg("stop") = static_cast<real_t>(5.0), py::arg("num_basis") = 8)
      .def_readwrite("start", &GaussianRBFConfig::start, "Start center distance in Angstroms.")
      .def_readwrite("stop", &GaussianRBFConfig::stop, "Stop center distance in Angstroms.")
      .def_readwrite("num_basis", &GaussianRBFConfig::num_basis, "Number of Gaussian basis centers.");

  // ------------------------------------------------------------------
  // PeriodicGraphBuilder
  // ------------------------------------------------------------------
  py::class_<PeriodicGraphBuilder>(mod, "PeriodicGraphBuilder",
                                   "Constructs periodic neighbor graphs for atomic GNN evaluation.")
      .def_static("build_graph", &PeriodicGraphBuilder::buildGraph, py::arg("cell"),
                  py::arg("cutoff_radius") = static_cast<real_t>(5.0), py::arg("include_self_loops") = false,
                  "Build periodic neighbor graph data for a unit cell.")
      .def_static("get_atomic_number", &PeriodicGraphBuilder::getAtomicNumber, py::arg("symbol"),
                  "Return atomic number (Z) for element symbol.")
      .def_static("compute_cutoff_envelope", &PeriodicGraphBuilder::computeCutoffEnvelope, py::arg("distance"),
                  py::arg("cutoff_radius"), "Compute smooth polynomial cutoff envelope.")
      .def_static("compute_bessel_basis", &PeriodicGraphBuilder::computeBesselBasis, py::arg("distance"),
                  py::arg("cutoff_radius"), py::arg("num_basis"), "Compute spherical Bessel radial basis.")
      .def_static(
          "compute_gaussian_rbf",
          [](real_t distance, real_t start, real_t stop, size_t num_basis) {
            return PeriodicGraphBuilder::computeGaussianRBF(distance,
                                                            {.start = start, .stop = stop, .num_basis = num_basis});
          },
          py::arg("distance"), py::arg("start"), py::arg("stop"), py::arg("num_basis"),
          "Compute Gaussian radial basis functions.")
      .def_static(
          "compute_gaussian_rbf",
          [](real_t distance, const GaussianRBFConfig &config) {
            return PeriodicGraphBuilder::computeGaussianRBF(distance, config);
          },
          py::arg("distance"), py::arg("config"),
          "Compute Gaussian radial basis functions using configuration struct.");

  // Free function alias for convenient top-level usage
  mod.def("build_periodic_graph", &PeriodicGraphBuilder::buildGraph, py::arg("cell"),
          py::arg("cutoff_radius") = static_cast<real_t>(5.0), py::arg("include_self_loops") = false,
          "Convenience helper to construct periodic neighbor graph for GNN evaluation.");

  // ------------------------------------------------------------------
  // MLIPInterface
  // ------------------------------------------------------------------
  py::class_<MLIPInterface, std::unique_ptr<MLIPInterface, py::nodelete>>(mod, "MLIPInterface",
                                                                          "Abstract interface for MLIP engines.")
      .def("get_model_name", &MLIPInterface::getModelName, "Return model descriptor name.")
      .def("evaluate", &MLIPInterface::evaluate, py::arg("cell"), "Evaluate model on an atomic cell.");

  // ------------------------------------------------------------------
  // MLIPOutput
  // ------------------------------------------------------------------
  py::class_<MLIPOutput>(mod, "MLIPOutput", "Container for machine learning interatomic potential outputs.")
      .def(py::init<>())
      .def_readwrite("total_energy", &MLIPOutput::total_energy, "Total predicted potential energy.")
      .def_readwrite("per_atom_energy", &MLIPOutput::per_atom_energy, "Site-resolved per-atom energy.")
      .def_readwrite("ldos", &MLIPOutput::ldos, "Local Density of States matrix [N_atoms x N_bins].")
      .def_readwrite("ldos_bins", &MLIPOutput::ldos_bins, "Number of LDOS energy bins.")
      .def_property_readonly(
          "forces",
          [](const MLIPOutput &out) -> py::array_t<real_t> {
            const size_t n_atoms = out.forces.size();
            py::array_t<real_t> arr({static_cast<py::ssize_t>(n_atoms), py::ssize_t(3)});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n_atoms; ++i) {
              buf(i, 0) = out.forces[i].x();
              buf(i, 1) = out.forces[i].y();
              buf(i, 2) = out.forces[i].z();
            }
            return arr;
          },
          "Predicted atomic forces as a (N, 3) NumPy array.");

  // ------------------------------------------------------------------
  // MLIPCalculator
  // ------------------------------------------------------------------
  py::class_<MLIPCalculator, BaseCalculator>(mod, "MLIPCalculator",
                                             "Calculator interface for machine learning interatomic potentials.")
      .def(py::init<>())
      .def_static(
          "calculate",
          [](const correlation::core::Cell &cell, const MLIPInterface *model) {
            return MLIPCalculator::calculate(cell, model);
          },
          py::arg("cell"), py::arg("model") = nullptr, "Evaluate MLIP on a given cell.");
}
