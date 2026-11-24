/**
 * @file mlip_bindings.cpp
 * @brief Pybind11 bindings for machine learning interatomic potentials (MLIP) and GNN graph
 * builder.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "mlip_bindings.hpp"
#include "calculators/ml/ChargeCorrelationCalculator.hpp"
#include "calculators/ml/EliashbergCalculator.hpp"
#include "calculators/ml/MLIPCalculator.hpp"
#include "calculators/ml/SpinCorrelationCalculator.hpp"
#include "calculators/scattering/LDOSPostProcessor.hpp"
#include "calculators/scattering/TDOSCalculator.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "mlip/GraphDescriptors.hpp"
#include "mlip/MLIPInterface.hpp"
#include "mlip/MaceGNNModel.hpp"
#include "mlip/OrbV3GNNModel.hpp"
#include "mlip/PeriodicGraphBuilder.hpp"

#include <span>

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

template <typename T> [[nodiscard]] py::array_t<T> makeOwned1D(std::span<const T> buffer) {
  return py::array_t<T>(static_cast<py::ssize_t>(buffer.size()), buffer.data());
}

void bindPeriodicGraphData(const py::module_ &mod) {
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
      .def(py::init<real_t, size_t, size_t, size_t, bool, bool>(),
           py::arg("r_max") = static_cast<real_t>(6.0), py::arg("num_rbf") = static_cast<size_t>(8),
           py::arg("l_max") = static_cast<size_t>(3),
           py::arg("max_num_neighbors") = static_cast<size_t>(120),
           py::arg("include_self_loops") = false, py::arg("compute_orb_features") = true)
      .def_readwrite("r_max", &OrbDescriptorConfig::r_max, "Radial cutoff distance in Angstroms.")
      .def_readwrite("num_rbf", &OrbDescriptorConfig::num_rbf,
                     "Number of Bessel radial basis functions.")
      .def_readwrite("l_max", &OrbDescriptorConfig::l_max, "Maximum spherical harmonics degree l.")
      .def_readwrite("max_num_neighbors", &OrbDescriptorConfig::max_num_neighbors,
                     "Maximum nearest periodic neighbors per atom.")
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
                  py::arg("max_num_neighbors") = static_cast<size_t>(0),
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
          py::arg("max_num_neighbors") = static_cast<size_t>(0),
          "Convenience helper to construct periodic neighbor graph for GNN evaluation.");
  mod.def(
      "build_orb_graph", &PeriodicGraphBuilder::buildOrbGraph, py::arg("cell"),
      py::arg("config") = OrbDescriptorConfig{},
      "Convenience helper to construct periodic neighbor graph and extract ORB-v3 descriptors.");
}

class PyMLIPInterface : public MLIPInterface {
public:
  using MLIPInterface::MLIPInterface;

  [[nodiscard]] std::string getModelName() const override {
    PYBIND11_OVERRIDE_PURE_NAME(std::string, MLIPInterface, "get_model_name", getModelName);
  }

  [[nodiscard]] MLIPOutput evaluate(const correlation::core::Cell &cell) const override {
    PYBIND11_OVERRIDE_PURE(MLIPOutput, MLIPInterface, evaluate, cell);
  }
};

[[nodiscard]] py::array_t<real_t>
vector3ListToNumpy(std::span<const correlation::math::Vector3<real_t>> vecs) {
  const size_t count = vecs.size();
  py::array_t<real_t> arr({static_cast<py::ssize_t>(count), static_cast<py::ssize_t>(3)});
  auto buf = arr.mutable_unchecked<2>();
  for (size_t idx = 0; idx < count; ++idx) {
    const auto row = static_cast<py::ssize_t>(idx);
    buf(row, 0) = vecs[idx].x();
    buf(row, 1) = vecs[idx].y();
    buf(row, 2) = vecs[idx].z();
  }
  return arr;
}

void numpyToVector3List(const py::array_t<real_t> &arr,
                        std::vector<correlation::math::Vector3<real_t>> &out) {
  const auto buf = arr.unchecked<2>();
  const auto count = static_cast<size_t>(buf.shape(0));
  out.resize(count);
  for (size_t idx = 0; idx < count; ++idx) {
    const auto row = static_cast<py::ssize_t>(idx);
    out[idx] = correlation::math::Vector3<real_t>{buf(row, 0), buf(row, 1), buf(row, 2)};
  }
}

[[nodiscard]] py::array_t<real_t> getLdosNumpy(const MLIPOutput &out) {
  if (out.ldos_bins == 0 || out.ldos.empty()) {
    return {};
  }
  const size_t n_atoms = out.ldos.size() / out.ldos_bins;
  py::array_t<real_t> arr(
      {static_cast<py::ssize_t>(n_atoms), static_cast<py::ssize_t>(out.ldos_bins)});
  auto buf = arr.mutable_unchecked<2>();
  for (size_t i = 0; i < n_atoms; ++i) {
    for (size_t bin_idx = 0; bin_idx < out.ldos_bins; ++bin_idx) {
      buf(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(bin_idx)) =
          out.ldos[i * out.ldos_bins + bin_idx];
    }
  }
  return arr;
}

void setLdosNumpy(MLIPOutput &out, const py::array_t<real_t> &arr) {
  const auto buf = arr.unchecked<2>();
  const auto n_atoms = static_cast<size_t>(buf.shape(0));
  const auto n_bins = static_cast<size_t>(buf.shape(1));
  out.ldos_bins = n_bins;
  out.ldos.resize(n_atoms * n_bins);
  for (size_t i = 0; i < n_atoms; ++i) {
    for (size_t bin_idx = 0; bin_idx < n_bins; ++bin_idx) {
      out.ldos[i * n_bins + bin_idx] =
          buf(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(bin_idx));
    }
  }
}

[[nodiscard]] py::array_t<real_t> getStressNumpy(const MLIPOutput &out) {
  py::array_t<real_t> arr({3, 3});
  auto buf = arr.mutable_unchecked<2>();
  for (size_t row = 0; row < 3; ++row) {
    for (size_t col = 0; col < 3; ++col) {
      buf(static_cast<py::ssize_t>(row), static_cast<py::ssize_t>(col)) = out.stress(row, col);
    }
  }
  return arr;
}

void setStressNumpy(MLIPOutput &out, const py::array_t<real_t> &arr) {
  const auto buf = arr.unchecked<2>();
  for (size_t row = 0; row < 3; ++row) {
    for (size_t col = 0; col < 3; ++col) {
      out.stress(row, col) = buf(static_cast<py::ssize_t>(row), static_cast<py::ssize_t>(col));
    }
  }
}

void bindMlipOutput(const py::module_ &mod) {
  py::class_<MLIPOutput>(mod, "MLIPOutput",
                         "Container for machine learning interatomic potential outputs.")
      .def(py::init<>())
      .def_readwrite("total_energy", &MLIPOutput::total_energy, "Total predicted potential energy.")
      .def_readwrite("per_atom_energy", &MLIPOutput::per_atom_energy,
                     "Site-resolved per-atom energy.")
      .def_property("ldos", &getLdosNumpy, &setLdosNumpy,
                    "Local Density of States matrix [N_atoms x N_bins].")
      .def_readwrite("ldos_bins", &MLIPOutput::ldos_bins, "Number of LDOS energy bins.")
      .def_readwrite("ldos_e_min", &MLIPOutput::ldos_e_min,
                     "Lower bound of LDOS energy grid in eV.")
      .def_readwrite("ldos_e_max", &MLIPOutput::ldos_e_max,
                     "Upper bound of LDOS energy grid in eV.")
      .def_readwrite("ldos_sigma", &MLIPOutput::ldos_sigma, "Gaussian broadening sigma in eV.")
      .def_property(
          "forces", [](const MLIPOutput &out) { return vector3ListToNumpy(out.forces); },
          [](MLIPOutput &out, const py::array_t<real_t> &arr) {
            numpyToVector3List(arr, out.forces);
          },
          "Predicted atomic forces as a (N, 3) NumPy array.")
      .def_property("stress", &getStressNumpy, &setStressNumpy,
                    "Virial stress tensor as a (3, 3) NumPy array.")
      .def_property(
          "charges",
          [](const MLIPOutput &out) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(out.charges);
          },
          [](MLIPOutput &out, const py::array_t<real_t> &arr) {
            const auto buf = arr.unchecked<1>();
            out.charges.resize(static_cast<size_t>(buf.shape(0)));
            for (size_t i = 0; i < out.charges.size(); ++i) {
              out.charges[i] = buf(static_cast<py::ssize_t>(i));
            }
          },
          "Per-atom partial or Born effective charges (e).")
      .def_property(
          "magmoms", [](const MLIPOutput &out) { return vector3ListToNumpy(out.magmoms); },
          [](MLIPOutput &out, const py::array_t<real_t> &arr) {
            numpyToVector3List(arr, out.magmoms);
          },
          "Per-atom site magnetic moment vectors (mu_B).")
      .def_readwrite("energy_variance", &MLIPOutput::energy_variance,
                     "Committee or ensemble energy prediction variance (eV^2).")
      .def_readwrite("force_variances", &MLIPOutput::force_variances,
                     "Committee or ensemble per-atom force variances ((eV/A)^2).");
}

void bindSpinChargeCalculators(const py::module_ &mod) {
  py::class_<SpinCorrelationResult>(
      mod, "SpinCorrelationResult",
      "Result of spatial spin-spin and magnetic moment correlation analysis.")
      .def(py::init<>())
      .def_property_readonly(
          "r_bins",
          [](const SpinCorrelationResult &res) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(res.r_bins);
          },
          "Radial distance bin centers in Angstroms.")
      .def_property_readonly(
          "spin_correlation",
          [](const SpinCorrelationResult &res) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(res.spin_correlation);
          },
          "Normalized spatial spin-spin alignment <S_i · S_j>(r) in [-1, 1].")
      .def_property_readonly(
          "magnetic_moment_corr",
          [](const SpinCorrelationResult &res) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(res.magnetic_moment_corr);
          },
          "Unnormalized magnetic moment correlation <m_i · m_j>(r) in mu_B^2.")
      .def_property_readonly(
          "pair_counts",
          [](const SpinCorrelationResult &res) -> py::array_t<uint64_t> {
            return makeOwned1D<uint64_t>(res.pair_counts);
          },
          "Pair observation counts per radial bin.");

  py::class_<SpinCorrelationCalculator, BaseCalculator>(
      mod, "SpinCorrelationCalculator",
      "Calculates spatial spin-spin correlation <S_i · S_j>(r) and magnetic RDF.")
      .def(py::init<>())
      .def_static(
          "calculate",
          [](const correlation::core::Cell &cell, const py::array_t<real_t> &magmoms_arr,
             real_t r_max, real_t r_bin_width) {
            const auto buf = magmoms_arr.unchecked<2>();
            const auto count = static_cast<size_t>(buf.shape(0));
            std::vector<correlation::math::Vector3<real_t>> vecs(count);
            for (size_t i = 0; i < count; ++i) {
              const auto row = static_cast<py::ssize_t>(i);
              vecs[i] = correlation::math::Vector3<real_t>{buf(row, 0), buf(row, 1), buf(row, 2)};
            }
            return SpinCorrelationCalculator::calculate(cell, vecs, r_max, r_bin_width);
          },
          py::arg("cell"), py::arg("magmoms"), py::arg("r_max") = static_cast<real_t>(10.0),
          py::arg("r_bin_width") = static_cast<real_t>(0.05),
          "Compute spatial spin correlation from Cell and (N, 3) magnetic moments array.")
      .def_static("calculate_from_mlip", &SpinCorrelationCalculator::calculateFromMLIP,
                  py::arg("cell"), py::arg("mlip_out"),
                  py::arg("r_max") = static_cast<real_t>(10.0),
                  py::arg("r_bin_width") = static_cast<real_t>(0.05),
                  "Compute spatial spin correlation from Cell and MLIPOutput.");

  py::class_<ChargeCorrelationResult>(
      mod, "ChargeCorrelationResult",
      "Result of charge-weighted pair distribution and screening analysis.")
      .def(py::init<>())
      .def_property_readonly(
          "r_bins",
          [](const ChargeCorrelationResult &res) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(res.r_bins);
          },
          "Radial distance bin centers in Angstroms.")
      .def_property_readonly(
          "g_qq",
          [](const ChargeCorrelationResult &res) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(res.g_qq);
          },
          "Charge-weighted pair distribution function g_qq(r).")
      .def_property_readonly(
          "charge_integral",
          [](const ChargeCorrelationResult &res) -> py::array_t<real_t> {
            return makeOwned1D<real_t>(res.charge_integral);
          },
          "Running cumulative charge neutrality screening integral Q(r).")
      .def_property_readonly(
          "pair_counts",
          [](const ChargeCorrelationResult &res) -> py::array_t<uint64_t> {
            return makeOwned1D<uint64_t>(res.pair_counts);
          },
          "Pair observation counts per radial bin.");

  py::class_<ChargeCorrelationCalculator, BaseCalculator>(
      mod, "ChargeCorrelationCalculator",
      "Calculates charge-weighted pair distribution g_qq(r) and screening integral Q(r).")
      .def(py::init<>())
      .def_static(
          "calculate",
          [](const correlation::core::Cell &cell, const py::array_t<real_t> &charges_arr,
             real_t r_max, real_t r_bin_width) {
            const auto buf = charges_arr.unchecked<1>();
            const auto count = static_cast<size_t>(buf.shape(0));
            std::vector<real_t> charges(count);
            for (size_t i = 0; i < count; ++i) {
              charges[i] = buf(static_cast<py::ssize_t>(i));
            }
            return ChargeCorrelationCalculator::calculate(cell, charges, r_max, r_bin_width);
          },
          py::arg("cell"), py::arg("charges"), py::arg("r_max") = static_cast<real_t>(10.0),
          py::arg("r_bin_width") = static_cast<real_t>(0.05),
          "Compute charge pair distribution g_qq(r) from Cell and 1D charges array.")
      .def_static("calculate_from_mlip", &ChargeCorrelationCalculator::calculateFromMLIP,
                  py::arg("cell"), py::arg("mlip_out"),
                  py::arg("r_max") = static_cast<real_t>(10.0),
                  py::arg("r_bin_width") = static_cast<real_t>(0.05),
                  "Compute charge pair distribution g_qq(r) from Cell and MLIPOutput.");
}

void bindMlipInterface(const py::module_ &mod) {
  py::class_<MLIPInterface, PyMLIPInterface>(mod, "MLIPInterface",
                                             "Abstract interface for MLIP engines.")
      .def(py::init<>())
      .def("get_model_name", &MLIPInterface::getModelName, "Return model descriptor name.")
      .def("evaluate", &MLIPInterface::evaluate, py::arg("cell"),
           "Evaluate model on an atomic cell.");

  bindMlipOutput(mod);

  py::class_<MaceGNNModel, MLIPInterface>(
      mod, "MaceGNNModel",
      "Native LibTorch MACE (Higher-Order Equivariant Message Passing) inference engine.")
      .def(py::init<std::string, std::string, real_t>(), py::arg("model_path"),
           py::arg("device") = "cpu", py::arg("cutoff") = static_cast<real_t>(5.0),
           "Construct MaceGNNModel from TorchScript .pt or .model file.")
      .def("is_loaded", &MaceGNNModel::isLoaded, "Check if model weights are loaded.")
      .def_property("cutoff", &MaceGNNModel::cutoff, &MaceGNNModel::setCutoff,
                    "Neighbor cutoff radius in Angstroms.");

  py::class_<OrbV3GNNModel, MLIPInterface>(
      mod, "OrbV3GNNModel", "Native C++23 ORB-v3 Graph Neural Network + LDOSHead inference engine.")
      .def(py::init<std::string>(), py::arg("model_path") = "",
           "Construct OrbV3GNNModel from binary weight file (auto-discovers default if empty).")
      .def_static("find_default_model_path", &OrbV3GNNModel::findDefaultModelPath,
                  py::arg("preferred_path") = "", "Find default cellulose_ldos.bin model path.")
      .def("is_loaded", &OrbV3GNNModel::isLoaded, "Check if model weights are loaded.")
      .def("model_path", &OrbV3GNNModel::modelPath, "Return resolved path to loaded weight file.");

  py::class_<MLIPCalculator, BaseCalculator>(
      mod, "MLIPCalculator", "Calculator interface for machine learning interatomic potentials.")
      .def(py::init<>())
      .def_static(
          "calculate",
          [](const correlation::core::Cell &cell, const MLIPInterface *model) {
            return MLIPCalculator::calculate(cell, model);
          },
          py::arg("cell"), py::arg("model") = nullptr, "Evaluate MLIP on a given cell.");

  bindSpinChargeCalculators(mod);
}

void bindTdos(const py::module_ &mod) {
  py::class_<BandGapInfo>(mod, "BandGapInfo",
                          "Electronic band-edge and gap summary around E_F = 0 eV.")
      .def(py::init<>())
      .def_readwrite("band_gap", &BandGapInfo::band_gap)
      .def_readwrite("vbm_energy", &BandGapInfo::vbm_energy)
      .def_readwrite("cbm_energy", &BandGapInfo::cbm_energy)
      .def_readwrite("is_metallic", &BandGapInfo::is_metallic);

  py::class_<LDOSPostProcessor>(mod, "LDOSPostProcessor",
                                "Stateless electronic LDOS/TDOS post-processing utilities.")
      .def_static(
          "broaden_spectrum",
          [](const std::vector<real_t> &energies, const std::vector<real_t> &dos,
             correlation::analysis::LDOSBroadeningKernel kernel, real_t sigma) {
            return LDOSPostProcessor::broadenSpectrum(energies, dos, kernel, sigma);
          },
          py::arg("energies"), py::arg("dos"), py::arg("kernel"), py::arg("sigma"))
      .def_static(
          "rebin_spectrum",
          [](const std::vector<real_t> &src_e, const std::vector<real_t> &src_dos,
             const std::vector<real_t> &dst_e) {
            return LDOSPostProcessor::rebinSpectrum(src_e, src_dos, dst_e);
          },
          py::arg("src_energies"), py::arg("src_dos"), py::arg("dst_energies"))
      .def_static(
          "evaluate_fermi_dos",
          [](const std::vector<real_t> &energies, const std::vector<real_t> &dos,
             real_t electronic_temp_k) {
            return LDOSPostProcessor::evaluateFermiDos(energies, dos, electronic_temp_k);
          },
          py::arg("energies"), py::arg("dos"),
          py::arg("electronic_temp_k") = static_cast<real_t>(0.0))
      .def_static(
          "detect_band_gap",
          [](const std::vector<real_t> &energies, const std::vector<real_t> &dos,
             real_t threshold) {
            return LDOSPostProcessor::detectBandGap(energies, dos, threshold);
          },
          py::arg("energies"), py::arg("dos"), py::arg("threshold") = static_cast<real_t>(1e-3))
      .def_static(
          "compute_integrated_dos",
          [](const std::vector<real_t> &energies, const std::vector<real_t> &dos) {
            return LDOSPostProcessor::computeIntegratedDos(energies, dos);
          },
          py::arg("energies"), py::arg("dos"));

  // ------------------------------------------------------------------
  // TDOSParams & LDOSResult
  // ------------------------------------------------------------------
  py::class_<TDOSParams>(mod, "TDOSParams",
                         "Configuration parameters for Total Density of States (TDOS) calculation.")
      .def(py::init<>())
      .def(py::init([](real_t e_min, real_t e_max, const MLIPInterface *model) {
             TDOSParams params;
             params.e_min = e_min;
             params.e_max = e_max;
             params.model = model;
             return params;
           }),
           py::arg("e_min") = static_cast<real_t>(-30.0),
           py::arg("e_max") = static_cast<real_t>(15.0), py::arg("model") = nullptr)
      .def_readwrite("e_min", &TDOSParams::e_min, "Lower energy bound in eV (relative to E_F).")
      .def_readwrite("e_max", &TDOSParams::e_max, "Upper energy bound in eV (relative to E_F).")
      .def_readwrite("model", &TDOSParams::model, "MLIPInterface model pointer.")
      .def_readwrite("broadening_kernel", &TDOSParams::broadening_kernel)
      .def_readwrite("broadening_sigma", &TDOSParams::broadening_sigma)
      .def_readwrite("electronic_temp_k", &TDOSParams::electronic_temp_k)
      .def_readwrite("gap_threshold", &TDOSParams::gap_threshold)
      .def_readwrite("rebin_num_bins", &TDOSParams::rebin_num_bins)
      .def_readwrite("compute_idos", &TDOSParams::compute_idos)
      .def_readwrite("retain_per_atom_ldos", &TDOSParams::retain_per_atom_ldos);

  py::class_<LDOSResult>(mod, "LDOSResult",
                         "Comprehensive TDOS, IDOS, electronic diagnostics, and per-atom LDOS.")
      .def(py::init<>())
      .def_readonly("tdos", &LDOSResult::tdos)
      .def_readonly("idos", &LDOSResult::idos)
      .def_readonly("diagnostics", &LDOSResult::diagnostics)
      .def_readonly("ldos_bins", &LDOSResult::ldos_bins)
      .def_property_readonly(
          "per_atom_ldos",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &res = obj.cast<const LDOSResult &>();
            const size_t n_atoms =
                (res.ldos_bins > 0) ? (res.per_atom_ldos.size() / res.ldos_bins) : 0;
            return makeEmptyOr2D<real_t>(res.per_atom_ldos, n_atoms, res.ldos_bins, obj);
          },
          "Zero-copy (N_atoms, N_bins) per-atom LDOS NumPy array.");

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
            const py::gil_scoped_release release;
            return TDOSCalculator::calculate(cell, params);
          },
          py::arg("cell"), py::arg("params") = TDOSParams{},
          "Calculate Total Density of States for a single cell.")
      .def_static(
          "calculate_detailed",
          [](const correlation::core::Cell &cell, const TDOSParams &params) {
            const py::gil_scoped_release release;
            return TDOSCalculator::calculateDetailed(cell, params);
          },
          py::arg("cell"), py::arg("params") = TDOSParams{},
          "Calculate TDOS, IDOS, electronic diagnostics, and optional per-atom LDOS for a cell.")
      .def_static(
          "calculate_trajectory",
          [](const correlation::core::Trajectory &traj, const TDOSParams &params) {
            const py::gil_scoped_release release;
            return TDOSCalculator::calculateTrajectory(traj, params, nullptr);
          },
          py::arg("traj"), py::arg("params") = TDOSParams{},
          "Calculate frame-averaged Total Density of States for a trajectory.")
      .def_static(
          "calculate_trajectory_detailed",
          [](const correlation::core::Trajectory &traj, const TDOSParams &params) {
            const py::gil_scoped_release release;
            return TDOSCalculator::calculateTrajectoryDetailed(traj, params, nullptr);
          },
          py::arg("traj"), py::arg("params") = TDOSParams{},
          "Calculate trajectory-averaged TDOS, IDOS, and electronic diagnostics.");
}

void bindEliashberg(const py::module_ &mod) {
  py::class_<EliashbergParams>(mod, "EliashbergParams",
                               "Parameters for Hybrid MLIP-VDOS Eliashberg alpha^2F(omega).")
      .def(py::init<>())
      .def_readwrite("model", &EliashbergParams::model)
      .def_readwrite("mu_star", &EliashbergParams::mu_star)
      .def_readwrite("coupling_mode", &EliashbergParams::coupling_mode)
      .def_readwrite("displacement_step", &EliashbergParams::displacement_step)
      .def_readwrite("max_coupling_frames", &EliashbergParams::max_coupling_frames)
      .def_readwrite("samples_per_species", &EliashbergParams::samples_per_species)
      .def_readwrite("electronic_temp_k", &EliashbergParams::electronic_temp_k)
      .def_readwrite("omega_min_mev", &EliashbergParams::omega_min_mev);

  py::class_<SpeciesCouplingWeights>(mod, "SpeciesCouplingWeights",
                                     "Species-resolved MLIP deformation sensitivity proxy.")
      .def(py::init<>())
      .def_readonly("eta_by_species", &SpeciesCouplingWeights::eta_by_species)
      .def_readonly("fermi_dos_by_species", &SpeciesCouplingWeights::fermi_dos_by_species)
      .def_readonly("total_fermi_dos", &SpeciesCouplingWeights::total_fermi_dos);

  py::class_<EliashbergResult>(mod, "EliashbergResult",
                               "Output bundle for Eliashberg alpha^2F(omega) estimation.")
      .def(py::init<>())
      .def_readonly("spectral_function", &EliashbergResult::spectral_function)
      .def_readonly("diagnostics", &EliashbergResult::diagnostics)
      .def_readonly("coupling_weights", &EliashbergResult::coupling_weights);

  py::class_<EliashbergCalculator, BaseCalculator>(
      mod, "EliashbergCalculator",
      "Hybrid MLIP-VDOS Eliashberg spectral function alpha^2F(omega) and Allen-Dynes T_c "
      "estimator.")
      .def(py::init<>())
      .def(py::init<const MLIPInterface *>(), py::arg("model") = nullptr)
      .def("set_model", &EliashbergCalculator::setModel, py::arg("model"))
      .def("get_model", &EliashbergCalculator::getModel)
      .def_static(
          "compute_from_spectrum",
          [](const std::vector<real_t> &omega_mev, const std::vector<real_t> &a2f, real_t mu_star,
             real_t omega_min_mev) {
            return EliashbergCalculator::computeFromSpectrum(omega_mev, a2f, mu_star,
                                                             omega_min_mev);
          },
          py::arg("omega_mev"), py::arg("a2f"), py::arg("mu_star") = static_cast<real_t>(0.10),
          py::arg("omega_min_mev") = static_cast<real_t>(0.5),
          "Compute lambda, omega_log, omega_2, McMillan and Allen-Dynes T_c from alpha^2F(omega).")
      .def_static(
          "calculate",
          [](const correlation::core::Trajectory &traj,
             const correlation::analysis::Histogram &vdos_vib, const EliashbergParams &params) {
            const py::gil_scoped_release release;
            return EliashbergCalculator::calculate(traj, vdos_vib, params, nullptr, nullptr);
          },
          py::arg("traj"), py::arg("vdos_vib"), py::arg("params") = EliashbergParams{},
          "Estimate Eliashberg alpha^2F(omega) and superconducting diagnostics from trajectory and "
          "vibrational VDOS.");
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
  bindEliashberg(mod);
  bindGraphDescriptors(mod);
}
