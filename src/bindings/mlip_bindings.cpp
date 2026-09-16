/**
 * @file mlip_bindings.cpp
 * @brief Pybind11 bindings for machine learning interatomic potentials (MLIP) and GNN graph
 * builder.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip_bindings.hpp"
#include "calculators/MLIPCalculator.hpp"
#include "calculators/TDOSCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "mlip/GraphDescriptors.hpp"
#include "mlip/MLIPInterface.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using correlation::real_t;
using namespace correlation::mlip;
using namespace correlation::calculators;

namespace {

template <typename T>
[[nodiscard]] py::array_t<T> makeEmptyOr1D(std::span<const T> buffer, size_t count,
                                           py::object &obj) {
  if (buffer.empty() || count == 0) {
    return {};
  }
  const auto count_val = static_cast<py::ssize_t>(count);
  const auto stride = static_cast<py::ssize_t>(sizeof(T));
  return py::array_t<T>({count_val}, {stride}, buffer.data(), obj);
}

template <typename T>
[[nodiscard]] py::array_t<T> makeEmptyOr2D(std::span<const T> buffer, size_t rows, size_t cols,
                                           py::object &obj) {
  if (buffer.empty() || rows == 0 || cols == 0) {
    return {};
  }
  const auto num_rows = static_cast<py::ssize_t>(rows);
  const auto num_cols = static_cast<py::ssize_t>(cols);
  const auto stride_row = static_cast<py::ssize_t>(num_cols * sizeof(T));
  const auto stride_col = static_cast<py::ssize_t>(sizeof(T));
  return py::array_t<T>({num_rows, num_cols}, {stride_row, stride_col}, buffer.data(), obj);
}

void bindPeriodicGraphData(py::module_ &mod) {
  py::class_<PeriodicGraphData>(
      mod, "PeriodicGraphData",
      "Container for periodic neighbor graph tensors ready for GNN model inference.")
      .def(py::init<>())
      .def_readonly("atom_count", &PeriodicGraphData::atom_count, "Total number of atoms N.")
      .def_readonly("edge_count", &PeriodicGraphData::edge_count,
                    "Total number of directed edges E.")
      .def_property_readonly(
          "positions",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr2D<real_t>(graph_data.positions_flat, graph_data.atom_count, 3, obj);
          },
          "Zero-copy access to atomic Cartesian coordinates as a (N, 3) NumPy array.")
      .def_property_readonly(
          "atomic_numbers",
          [](py::object &obj) -> py::array_t<int64_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr1D<int64_t>(graph_data.atomic_numbers, graph_data.atom_count, obj);
          },
          "Zero-copy access to atomic numbers (Z) as a (N,) NumPy array.")
      .def_property_readonly(
          "edge_index",
          [](py::object &obj) -> py::array_t<int64_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            if (graph_data.edge_index_flat.empty() || graph_data.edge_count == 0) {
              return {};
            }
            const py::ssize_t rows = 2;
            const auto cols = static_cast<py::ssize_t>(graph_data.edge_count);
            const auto stride_row =
                static_cast<py::ssize_t>(graph_data.edge_count * sizeof(int64_t));
            const auto stride_col = static_cast<py::ssize_t>(sizeof(int64_t));
            return py::array_t<int64_t>({rows, cols}, {stride_row, stride_col},
                                        graph_data.edge_index_flat.data(), obj);
          },
          "Zero-copy access to directed edge indices (COO format) as a (2, E) NumPy array.")
      .def_property_readonly(
          "edge_shifts",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr2D<real_t>(graph_data.edge_shifts_flat, graph_data.edge_count, 3,
                                         obj);
          },
          "Zero-copy access to periodic cell integer shift vectors as a (E, 3) NumPy array.")
      .def_property_readonly(
          "edge_vectors",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr2D<real_t>(graph_data.edge_vectors_flat, graph_data.edge_count, 3,
                                         obj);
          },
          "Zero-copy access to Cartesian displacement vectors r_ij as a (E, 3) NumPy array.")
      .def_property_readonly(
          "edge_distances",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr1D<real_t>(graph_data.edge_distances, graph_data.edge_count, obj);
          },
          "Zero-copy access to Euclidean edge distances ||r_ij|| as a (E,) NumPy array.")
      .def_property_readonly(
          "cell",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr2D<real_t>(graph_data.cell_flat, 3, 3, obj);
          },
          "Zero-copy access to lattice matrix as a (3, 3) NumPy array.")
      .def_property_readonly(
          "edge_spherical_harmonics",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            const size_t cols =
                (graph_data.edge_count > 0)
                    ? (graph_data.edge_spherical_harmonics_flat.size() / graph_data.edge_count)
                    : 0;
            return makeEmptyOr2D<real_t>(graph_data.edge_spherical_harmonics_flat,
                                         graph_data.edge_count, cols, obj);
          },
          "Zero-copy access to spherical harmonics features as a (E, (l_max+1)^2) NumPy array.")
      .def_property_readonly(
          "edge_unit_vectors",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr2D<real_t>(graph_data.edge_unit_vectors_flat, graph_data.edge_count,
                                         3, obj);
          },
          "Zero-copy access to normalized unit displacement vectors r_hat_ij as a (E, 3) NumPy array.")
      .def_property_readonly(
          "edge_radial_basis",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            const size_t cols =
                (graph_data.edge_count > 0)
                    ? (graph_data.edge_radial_basis_flat.size() / graph_data.edge_count)
                    : 0;
            return makeEmptyOr2D<real_t>(graph_data.edge_radial_basis_flat, graph_data.edge_count,
                                         cols, obj);
          },
          "Zero-copy access to edge radial basis features as a (E, num_rbf) NumPy array.")
      .def_property_readonly(
          "edge_cutoff_envelope",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr1D<real_t>(graph_data.edge_cutoff_envelope_flat,
                                         graph_data.edge_count, obj);
          },
          "Zero-copy access to polynomial cutoff envelope values f_c(d) as a (E,) NumPy array.")
      .def_property_readonly(
          "edge_orb_features",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            const size_t cols =
                (graph_data.edge_count > 0)
                    ? (graph_data.edge_orb_features_flat.size() / graph_data.edge_count)
                    : 0;
            return makeEmptyOr2D<real_t>(graph_data.edge_orb_features_flat, graph_data.edge_count,
                                         cols, obj);
          },
          "Zero-copy access to fused ORB-v3 edge features f_cut * (RBF (x) Y_lm) as a (E, num_rbf * num_sh) NumPy array.")
      .def_property_readonly(
          "cna_labels",
          [](py::object &obj) -> py::array_t<int> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr1D<int>(graph_data.cna_labels, graph_data.cna_labels.size(), obj);
          },
          "Zero-copy access to CNA classification labels as a (N,) NumPy array.")
      .def_property_readonly(
          "coordination_desc",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            return makeEmptyOr1D<real_t>(graph_data.coordination_desc,
                                         graph_data.coordination_desc.size(), obj);
          },
          "Zero-copy access to coordination numbers as a (N,) NumPy array.")
      .def_property_readonly(
          "ring_desc",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &graph_data = obj.cast<const PeriodicGraphData &>();
            const size_t cols = (graph_data.atom_count > 0)
                                    ? (graph_data.ring_desc.size() / graph_data.atom_count)
                                    : 0;
            return makeEmptyOr2D<real_t>(graph_data.ring_desc, graph_data.atom_count, cols, obj);
          },
          "Zero-copy access to per-atom ring descriptors as a (N, max_ring_size) NumPy array.");
}

void bindPeriodicGraphBuilder(py::module_ &mod) {
  // ------------------------------------------------------------------
  // GaussianRBFConfig
  // ------------------------------------------------------------------
  py::class_<GaussianRBFConfig>(
      mod, "GaussianRBFConfig",
      "Configuration parameters for Gaussian radial basis function expansion.")
      .def(py::init<real_t, real_t, size_t>(), py::arg("start") = static_cast<real_t>(0.0),
           py::arg("stop") = static_cast<real_t>(5.0), py::arg("num_basis") = 8)
      .def_readwrite("start", &GaussianRBFConfig::start, "Start center distance in Angstroms.")
      .def_readwrite("stop", &GaussianRBFConfig::stop, "Stop center distance in Angstroms.")
      .def_readwrite("num_basis", &GaussianRBFConfig::num_basis,
                     "Number of Gaussian basis centers.");

  // ------------------------------------------------------------------
  // OrbDescriptorConfig
  // ------------------------------------------------------------------
  py::class_<OrbDescriptorConfig>(
      mod, "OrbDescriptorConfig",
      "Configuration parameters for ORB-v3 graph and descriptor extraction.")
      .def(py::init<real_t, size_t, size_t, bool, bool>(),
           py::arg("r_max") = static_cast<real_t>(6.0), py::arg("num_rbf") = static_cast<size_t>(8),
           py::arg("l_max") = static_cast<size_t>(3), py::arg("include_self_loops") = false,
           py::arg("compute_orb_features") = true)
      .def_readwrite("r_max", &OrbDescriptorConfig::r_max, "Radial cutoff distance in Angstroms.")
      .def_readwrite("num_rbf", &OrbDescriptorConfig::num_rbf,
                     "Number of Bessel radial basis functions.")
      .def_readwrite("l_max", &OrbDescriptorConfig::l_max, "Maximum spherical harmonics degree l.")
      .def_readwrite("include_self_loops", &OrbDescriptorConfig::include_self_loops,
                     "Whether to include self loops.")
      .def_readwrite("compute_orb_features", &OrbDescriptorConfig::compute_orb_features,
                     "Whether to compute fused outer-product [E, 128] ORB edge features.");

  // ------------------------------------------------------------------
  // PeriodicGraphBuilder
  // ------------------------------------------------------------------
  py::class_<PeriodicGraphBuilder>(mod, "PeriodicGraphBuilder",
                                   "Constructs periodic neighbor graphs for atomic GNN evaluation.")
      .def_static("build_graph", &PeriodicGraphBuilder::buildGraph, py::arg("cell"),
                  py::arg("cutoff_radius") = static_cast<real_t>(5.0),
                  py::arg("include_self_loops") = false, py::arg("l_max") = static_cast<size_t>(0),
                  "Build periodic neighbor graph data for a unit cell.")
      .def_static("build_orb_graph", &PeriodicGraphBuilder::buildOrbGraph, py::arg("cell"),
                  py::arg("config") = OrbDescriptorConfig{},
                  "Construct periodic neighbor graph and compute native ORB-v3 descriptors.")
      .def_static("get_atomic_number", &PeriodicGraphBuilder::getAtomicNumber, py::arg("symbol"),
                  "Return atomic number (Z) for element symbol.")
      .def_static("compute_cutoff_envelope", &PeriodicGraphBuilder::computeCutoffEnvelope,
                  py::arg("distance"), py::arg("cutoff_radius"),
                  "Compute smooth polynomial cutoff envelope.")
      .def_static("compute_orb_cutoff_envelope", &PeriodicGraphBuilder::computeOrbCutoffEnvelope,
                  py::arg("distance"), py::arg("cutoff_radius"),
                  "Compute ORB-v3 polynomial cutoff envelope of order p=4.")
      .def_static("compute_bessel_basis", &PeriodicGraphBuilder::computeBesselBasis,
                  py::arg("distance"), py::arg("cutoff_radius"), py::arg("num_basis"),
                  "Compute spherical Bessel radial basis.")
      .def_static("compute_orb_bessel_basis", &PeriodicGraphBuilder::computeOrbBesselBasis,
                  py::arg("distance"), py::arg("config") = OrbDescriptorConfig{},
                  "Compute ORB-v3 Bessel radial basis functions.")
      .def_static(
          "compute_orb_spherical_harmonics",
          [](const std::vector<real_t> &unit_vec, size_t l_max) {
            if (unit_vec.size() < 3) {
              throw std::invalid_argument("unit_vec must have at least 3 elements");
            }
            return PeriodicGraphBuilder::computeOrbSphericalHarmonics(
                correlation::math::Vector3<real_t>{unit_vec[0], unit_vec[1], unit_vec[2]}, l_max);
          },
          py::arg("unit_vec"), py::arg("l_max") = static_cast<size_t>(3),
          "Compute ORB-v3 e3nn component-normalized spherical harmonics.")
      .def_static(
          "compute_gaussian_rbf",
          [](real_t distance, real_t start, real_t stop, size_t num_basis) {
            return PeriodicGraphBuilder::computeGaussianRBF(
                distance, {.start = start, .stop = stop, .num_basis = num_basis});
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

  // Free function aliases for convenient top-level usage
  mod.def("build_periodic_graph", &PeriodicGraphBuilder::buildGraph, py::arg("cell"),
          py::arg("cutoff_radius") = static_cast<real_t>(5.0),
          py::arg("include_self_loops") = false, py::arg("l_max") = static_cast<size_t>(0),
          "Convenience helper to construct periodic neighbor graph for GNN evaluation.");
  mod.def(
      "build_orb_graph", &PeriodicGraphBuilder::buildOrbGraph, py::arg("cell"),
      py::arg("config") = OrbDescriptorConfig{},
      "Convenience helper to construct periodic neighbor graph and extract ORB-v3 descriptors.");
}

void bindMlipInterface(py::module_ &mod) {
  // ------------------------------------------------------------------
  // MLIPInterface
  // ------------------------------------------------------------------
  py::class_<MLIPInterface, std::unique_ptr<MLIPInterface, py::nodelete>>(
      mod, "MLIPInterface", "Abstract interface for MLIP engines.")
      .def("get_model_name", &MLIPInterface::getModelName, "Return model descriptor name.")
      .def("evaluate", &MLIPInterface::evaluate, py::arg("cell"),
           "Evaluate model on an atomic cell.");

  // ------------------------------------------------------------------
  // MLIPOutput
  // ------------------------------------------------------------------
  py::class_<MLIPOutput>(mod, "MLIPOutput",
                         "Container for machine learning interatomic potential outputs.")
      .def(py::init<>())
      .def_readwrite("total_energy", &MLIPOutput::total_energy, "Total predicted potential energy.")
      .def_readwrite("per_atom_energy", &MLIPOutput::per_atom_energy,
                     "Site-resolved per-atom energy.")
      .def_readwrite("ldos", &MLIPOutput::ldos,
                     "Local Density of States matrix [N_atoms x N_bins].")
      .def_readwrite("ldos_bins", &MLIPOutput::ldos_bins, "Number of LDOS energy bins.")
      .def_property_readonly(
          "forces",
          [](const MLIPOutput &out) -> py::array_t<real_t> {
            const size_t n_atoms = out.forces.size();
            py::array_t<real_t> arr(
                {static_cast<py::ssize_t>(n_atoms), static_cast<py::ssize_t>(3)});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t idx = 0; idx < n_atoms; ++idx) {
              buf(idx, 0) = out.forces[idx].x();
              buf(idx, 1) = out.forces[idx].y();
              buf(idx, 2) = out.forces[idx].z();
            }
            return arr;
          },
          "Predicted atomic forces as a (N, 3) NumPy array.");

  // ------------------------------------------------------------------
  // MLIPCalculator
  // ------------------------------------------------------------------
  py::class_<MLIPCalculator, BaseCalculator>(
      mod, "MLIPCalculator", "Calculator interface for machine learning interatomic potentials.")
      .def(py::init<>())
      .def_static(
          "calculate",
          [](const correlation::core::Cell &cell, const MLIPInterface *model) {
            return MLIPCalculator::calculate(cell, model);
          },
          py::arg("cell"), py::arg("model") = nullptr, "Evaluate MLIP on a given cell.");
}

void bindTdos(py::module_ &mod) {
  // ------------------------------------------------------------------
  // TDOSParams
  // ------------------------------------------------------------------
  py::class_<TDOSParams>(mod, "TDOSParams",
                         "Configuration parameters for Total Density of States (TDOS) calculation.")
      .def(py::init<real_t, real_t, const MLIPInterface *>(),
           py::arg("e_min") = static_cast<real_t>(-15.0),
           py::arg("e_max") = static_cast<real_t>(5.0), py::arg("model") = nullptr)
      .def_readwrite("e_min", &TDOSParams::e_min, "Lower energy bound in eV (relative to E_F).")
      .def_readwrite("e_max", &TDOSParams::e_max, "Upper energy bound in eV (relative to E_F).")
      .def_readwrite("model", &TDOSParams::model, "MLIPInterface model pointer.");

  // ------------------------------------------------------------------
  // TDOSCalculator
  // ------------------------------------------------------------------
  py::class_<TDOSCalculator, BaseCalculator>(mod, "TDOSCalculator",
                                             "Calculator for MLIP Total Density of States (TDOS).")
      .def(py::init<>())
      .def(py::init<const MLIPInterface *>(), py::arg("model") = nullptr,
           "Construct TDOSCalculator with an optional MLIP model engine.")
      .def("set_model", &TDOSCalculator::setModel, py::arg("model"),
           "Set or update the MLIP model engine.")
      .def("get_model", &TDOSCalculator::getModel, "Get pointer to attached MLIP model engine.")
      .def_static(
          "calculate",
          [](const correlation::core::Cell &cell, const TDOSParams &params) {
            return TDOSCalculator::calculate(cell, params);
          },
          py::arg("cell"), py::arg("params") = TDOSParams{},
          "Calculate Total Density of States for a single cell.")
      .def_static(
          "calculate_trajectory",
          [](const correlation::core::Trajectory &traj, const TDOSParams &params) {
            return TDOSCalculator::calculateTrajectory(traj, params, nullptr);
          },
          py::arg("traj"), py::arg("params") = TDOSParams{},
          "Calculate frame-averaged Total Density of States for a trajectory.");
}

void bindGraphDescriptors(py::module_ &mod) {
  py::enum_<CNALabel>(mod, "CNALabel", "Canonical Common Neighbor Analysis structural motifs.")
      .value("OTHER", CNALabel::Other)
      .value("FCC", CNALabel::FCC)
      .value("HCP", CNALabel::HCP)
      .value("BCC", CNALabel::BCC)
      .value("ICO", CNALabel::ICO)
      .export_values();

  py::class_<GraphDescriptors>(
      mod, "GraphDescriptors",
      "Extracts topological, structural, and spectral descriptors from PeriodicGraphData.")
      .def_static("compute_ring_statistics_descriptor",
                  &GraphDescriptors::computeRingStatisticsDescriptor, py::arg("graph"),
                  py::arg("max_size") = 6,
                  "Compute per-atom ring statistics embedding using cycle basis detection.")
      .def_static("compute_cna_descriptor", &GraphDescriptors::computeCNADescriptor,
                  py::arg("graph"),
                  "Compute per-atom Common Neighbor Analysis (CNA) classification labels.")
      .def_static("compute_coordination_embedding", &GraphDescriptors::computeCoordinationEmbedding,
                  py::arg("graph"), "Compute per-atom coordination number embedding.")
      .def_static("compute_graph_spectrum", &GraphDescriptors::computeGraphSpectrum,
                  py::arg("graph"), py::arg("k"),
                  "Compute top-k eigenvalues of the graph adjacency matrix.")
      .def_static("populate_descriptors", &GraphDescriptors::populateDescriptors, py::arg("graph"),
                  py::arg("max_ring_size") = 6,
                  "Populate all descriptor fields in PeriodicGraphData in-place.");

  mod.def("compute_ring_statistics_descriptor", &GraphDescriptors::computeRingStatisticsDescriptor,
          py::arg("graph"), py::arg("max_size") = 6,
          "Compute per-atom ring statistics embedding using cycle basis detection.");
  mod.def("compute_cna_descriptor", &GraphDescriptors::computeCNADescriptor, py::arg("graph"),
          "Compute per-atom Common Neighbor Analysis (CNA) classification labels.");
  mod.def("compute_coordination_embedding", &GraphDescriptors::computeCoordinationEmbedding,
          py::arg("graph"), "Compute per-atom coordination number embedding.");
  mod.def("compute_graph_spectrum", &GraphDescriptors::computeGraphSpectrum, py::arg("graph"),
          py::arg("k"), "Compute top-k eigenvalues of the graph adjacency matrix.");
  mod.def("populate_descriptors", &GraphDescriptors::populateDescriptors, py::arg("graph"),
          py::arg("max_ring_size") = 6,
          "Populate all descriptor fields in PeriodicGraphData in-place.");
}

} // namespace

void init_mlip(py::module_ &mod) {
  bindPeriodicGraphData(mod);
  bindPeriodicGraphBuilder(mod);
  bindMlipInterface(mod);
  bindTdos(mod);
  bindGraphDescriptors(mod);
}
