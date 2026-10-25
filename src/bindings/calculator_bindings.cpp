/**
 * @file calculator_bindings.cpp
 * @brief Pybind11 bindings for correlation::calculators types.
 *
 * Exposes:
 *   - BaseCalculator   (abstract base with Python trampoline)
 *   - list_calculators()     (module-level free function)
 *   - get_calculator()       (module-level free function)
 *
 * Note: CalculatorFactory is a move-only singleton holding unique_ptrs.
 * It is not exposed as a Python class; module-level free functions provide
 * access to its registered calculators.
 */

#include "calculator_bindings.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "calculators/BaseCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/order/LocalEntropyCalculator.hpp"
#include "calculators/order/TenWoldeFrenkelCalculator.hpp"
#include "calculators/order/TetrahedralCalculator.hpp"
#include "calculators/spatial/G6Calculator.hpp"
#include "core/Cell.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <span>
#include <stdexcept>

namespace py = pybind11;
using correlation::real_t;
using namespace correlation::calculators;
using namespace correlation::analysis;
using namespace correlation::core;

// ------------------------------------------------------------------
// Trampoline — lets users subclass BaseCalculator in Python
// ------------------------------------------------------------------
namespace {
class PyBaseCalculator : public BaseCalculator {
  mutable std::string cached_name_;
  mutable std::string cached_short_name_;
  mutable std::string cached_group_;
  mutable std::string cached_description_;

public:
  using BaseCalculator::BaseCalculator;

  std::string_view getName() const override {
    const py::gil_scoped_acquire gil;
    const py::function overload = py::get_overload(this, "get_name");
    if (overload) {
      cached_name_ = overload().cast<std::string>();
      return cached_name_;
    }
    py::pybind11_fail("Tried to call pure virtual function \"BaseCalculator::getName\"");
  }

  std::string_view getShortName() const override {
    const py::gil_scoped_acquire gil;
    const py::function overload = py::get_overload(this, "get_short_name");
    if (overload) {
      cached_short_name_ = overload().cast<std::string>();
      return cached_short_name_;
    }
    py::pybind11_fail("Tried to call pure virtual function \"BaseCalculator::getShortName\"");
  }

  std::string_view getGroup() const override {
    const py::gil_scoped_acquire gil;
    const py::function overload = py::get_overload(this, "get_group");
    if (overload) {
      cached_group_ = overload().cast<std::string>();
      return cached_group_;
    }
    py::pybind11_fail("Tried to call pure virtual function \"BaseCalculator::getGroup\"");
  }

  std::string_view getDescription() const override {
    const py::gil_scoped_acquire gil;
    const py::function overload = py::get_overload(this, "get_description");
    if (overload) {
      cached_description_ = overload().cast<std::string>();
      return cached_description_;
    }
    py::pybind11_fail("Tried to call pure virtual function \"BaseCalculator::getDescription\"");
  }

  bool isFrameCalculator() const override {
    PYBIND11_OVERRIDE_PURE(bool, BaseCalculator, isFrameCalculator);
  }
  bool isTrajectoryCalculator() const override {
    PYBIND11_OVERRIDE_PURE(bool, BaseCalculator, isTrajectoryCalculator);
  }
  bool isConfigured() const override { PYBIND11_OVERRIDE(bool, BaseCalculator, isConfigured); }
  void calculateFrame(DistributionFunctions &dists,
                      const AnalysisSettings &settings) const override {
    PYBIND11_OVERRIDE(void, BaseCalculator, calculateFrame, dists, settings);
  }
  void calculateTrajectory(DistributionFunctions &dists, const correlation::core::Trajectory &traj,
                           const AnalysisSettings &settings) const override {
    PYBIND11_OVERRIDE(void, BaseCalculator, calculateTrajectory, dists, traj, settings);
  }
};

template <typename T> py::array_t<T> toNumpyVector(const std::vector<T> &vec) {
  py::array_t<T> arr(static_cast<py::ssize_t>(vec.size()));
  if (!vec.empty()) {
    std::copy(vec.begin(), vec.end(), arr.mutable_data());
  }
  return arr;
}

} // namespace

void init_calculators(py::module_ &mod) {
  // ------------------------------------------------------------------
  // BaseCalculator
  // ------------------------------------------------------------------
  py::class_<BaseCalculator, PyBaseCalculator>(
      mod, "BaseCalculator",
      "Abstract base class for all analysis calculators.\n\n"
      "Subclass this in Python to create custom calculators.")
      .def(py::init<>())
      .def("get_name", &BaseCalculator::getName,
           "Full display name of the calculator (e.g. 'g(r), J(r), G(r)').")
      .def("get_short_name", &BaseCalculator::getShortName,
           "Short identifier used as a key (e.g. 'RDF').")
      .def("get_group", &BaseCalculator::getGroup,
           "UI group this calculator belongs to (e.g. 'Radial', 'Angular').")
      .def("get_description", &BaseCalculator::getDescription,
           "Human-readable description of the calculator.")
      .def("is_frame_calculator", &BaseCalculator::isFrameCalculator,
           "True if this calculator operates per-frame.")
      .def("is_trajectory_calculator", &BaseCalculator::isTrajectoryCalculator,
           "True if this calculator operates over the full trajectory.")
      .def(
          "is_configured", &BaseCalculator::isConfigured,
          "True if the calculator has all required dependencies and models configured to execute.");

  // ------------------------------------------------------------------
  // CalculatorFactory access — module-level free functions
  // (CalculatorFactory holds unique_ptrs and cannot be copy-constructed)
  // ------------------------------------------------------------------
  mod.def(
      "list_calculators",
      [] {
        std::vector<std::string> names;
        for (const auto &calculator : CalculatorFactory::instance().getCalculators()) {
          names.emplace_back(calculator->getShortName());
        }
        return names;
      },
      "Return a list of short names for all registered calculators\n"
      "(e.g. ['RDF', 'SQ', 'PAD', ...]).");

  mod.def(
      "get_calculator",
      [](const std::string &name) -> const BaseCalculator * {
        const auto *calculator = CalculatorFactory::instance().getCalculator(name);
        if (!calculator) {
          throw std::runtime_error("No calculator registered with name: " + name);
        }
        return calculator;
      },
      py::arg("name"), py::return_value_policy::reference,
      "Look up a registered calculator by its short name (e.g. 'RDF').\n\n"
      "Returns a reference into the factory's singleton registry.\n\n"
      "Parameters\n----------\n"
      "name : str\n"
      "    Short name of the desired calculator.\n\n"
      "Raises\n------\n"
      "RuntimeError\n"
      "    If no calculator with that name is registered.");

  mod.def(
      "get_all_calculators",
      [] {
        std::vector<const BaseCalculator *> out;
        for (const auto &calculator : CalculatorFactory::instance().getCalculators()) {
          out.push_back(calculator.get());
        }
        return out;
      },
      py::return_value_policy::reference, "Return a list of all registered calculator objects.");

  // ------------------------------------------------------------------
  // Parameter types for Order & Disorder Calculators
  // ------------------------------------------------------------------
  py::class_<TenWoldeFrenkelParams>(mod, "TenWoldeFrenkelParams",
                                    "Parameters for Ten Wolde-Frenkel bond coherence analysis.")
      .def(py::init<real_t, size_t>(), py::arg("chi_cut") = 0.7, py::arg("nconn_min") = 7)
      .def_readwrite("chi_cut", &TenWoldeFrenkelParams::chi_cut,
                     "Coherence dot-product threshold (default 0.7).")
      .def_readwrite("nconn_min", &TenWoldeFrenkelParams::nconn_min,
                     "Minimum coherent neighbors to classify as solid-like (default 7).");

  py::class_<TetrahedralParams>(mod, "TetrahedralParams",
                                "Parameters for tetrahedral order parameter q_tet calculation.")
      .def(py::init<size_t, real_t, real_t>(), py::arg("num_bins") = 50,
           py::arg("min_domain") = 0.0, py::arg("max_domain") = 1.0)
      .def_readwrite("num_bins", &TetrahedralParams::num_bins, "Number of histogram bins.")
      .def_readwrite("min_domain", &TetrahedralParams::min_domain, "Lower domain bound.")
      .def_readwrite("max_domain", &TetrahedralParams::max_domain, "Upper domain bound.");

  py::class_<LocalEntropyParams>(mod, "LocalEntropyParams",
                                 "Parameters for Local Entropic Fingerprint (s2) calculation.")
      .def(py::init<real_t, real_t>(), py::arg("cutoff") = 5.0, py::arg("sigma") = 0.2)
      .def_readwrite("cutoff", &LocalEntropyParams::cutoff,
                     "Radial cutoff radius for entropy integration (default 5.0 A).")
      .def_readwrite("sigma", &LocalEntropyParams::sigma,
                     "Gaussian smoothing standard deviation (default 0.2 A).");

  py::class_<CorrelationLengthFit>(mod, "CorrelationLengthFit",
                                   "Fitted correlation length results from g6(r) profile.")
      .def(py::init<>())
      .def_readwrite("xi_6", &CorrelationLengthFit::xi_6,
                     "Correlation length in Angstroms (inf if crystalline plateau).")
      .def_readwrite("r_squared", &CorrelationLengthFit::r_squared,
                     "Coefficient of determination R^2 of log-linear regression fit.")
      .def_readwrite("is_plateau", &CorrelationLengthFit::is_plateau,
                     "True if flat crystalline plateau without exponential decay.");

  // ------------------------------------------------------------------
  // Free functions for Per-Atom & Structural Calculators
  // ------------------------------------------------------------------
  mod.def(
      "compute_ten_wolde_frenkel",
      [](const Cell &cell, const StructureAnalyzer *neighbors,
         const TenWoldeFrenkelParams &params) {
        const py::gil_scoped_release release;
        const auto res = TenWoldeFrenkelCalculator::calculatePerAtom(cell, neighbors, params);
        const py::gil_scoped_acquire acquire;
        return toNumpyVector(res);
      },
      py::arg("cell"), py::arg("neighbors") = nullptr, py::arg("params") = TenWoldeFrenkelParams(),
      "Computes per-atom Ten Wolde-Frenkel coherent neighbor count (N_conn).\n\n"
      "Parameters\n----------\n"
      "cell : Cell\n"
      "neighbors : StructureAnalyzer, optional\n"
      "params : TenWoldeFrenkelParams, optional\n\n"
      "Returns\n-------\n"
      "numpy.ndarray : 1D array of shape (N,) containing N_conn per atom.");

  mod.def(
      "compute_tetrahedral_order",
      [](const Cell &cell, const StructureAnalyzer *neighbors, const TetrahedralParams &params) {
        const py::gil_scoped_release release;
        const auto res = TetrahedralCalculator::calculatePerAtom(cell, neighbors, params);
        const py::gil_scoped_acquire acquire;
        return toNumpyVector(res);
      },
      py::arg("cell"), py::arg("neighbors") = nullptr, py::arg("params") = TetrahedralParams(),
      "Computes per-atom Chau-Hardwick tetrahedral order parameter (q_tet).\n\n"
      "Parameters\n----------\n"
      "cell : Cell\n"
      "neighbors : StructureAnalyzer, optional\n"
      "params : TetrahedralParams, optional\n\n"
      "Returns\n-------\n"
      "numpy.ndarray : 1D array of shape (N,) containing q_tet per atom.");

  mod.def(
      "compute_local_entropy",
      [](const Cell &cell, const StructureAnalyzer *neighbors, const LocalEntropyParams &params) {
        const py::gil_scoped_release release;
        const auto res = LocalEntropyCalculator::calculatePerAtom(cell, neighbors, params);
        const py::gil_scoped_acquire acquire;
        return toNumpyVector(res);
      },
      py::arg("cell"), py::arg("neighbors") = nullptr, py::arg("params") = LocalEntropyParams(),
      "Computes per-atom Local Entropic Fingerprint (s2).\n\n"
      "Parameters\n----------\n"
      "cell : Cell\n"
      "neighbors : StructureAnalyzer, optional\n"
      "params : LocalEntropyParams, optional\n\n"
      "Returns\n-------\n"
      "numpy.ndarray : 1D array of shape (N,) containing s2 per atom.");

  mod.def(
      "compute_g6",
      [](const Cell &cell, const StructureAnalyzer *neighbors, real_t r_max, real_t r_bin_width) {
        const py::gil_scoped_release release;
        auto hist = G6Calculator::calculate(cell, neighbors, r_max, r_bin_width);
        const py::gil_scoped_acquire acquire;
        return hist;
      },
      py::arg("cell"), py::arg("neighbors"), py::arg("r_max") = 10.0, py::arg("r_bin_width") = 0.05,
      "Computes radial bond-orientational correlation function g6(r).\n\n"
      "Parameters\n----------\n"
      "cell : Cell\n"
      "neighbors : StructureAnalyzer\n"
      "r_max : float, optional (default 10.0)\n"
      "r_bin_width : float, optional (default 0.05)\n\n"
      "Returns\n-------\n"
      "Histogram : 1D Histogram containing radial g6(r) profile.");

  mod.def(
      "fit_correlation_length",
      [](const py::array_t<real_t, py::array::c_style | py::array::forcecast> &r_bins,
         const py::array_t<real_t, py::array::c_style | py::array::forcecast> &g6_values) {
        if (r_bins.ndim() != 1 || g6_values.ndim() != 1) {
          throw std::invalid_argument("r_bins and g6_values must be 1D arrays.");
        }
        const std::span<const real_t> r_span(r_bins.data(), static_cast<size_t>(r_bins.size()));
        const std::span<const real_t> g6_span(g6_values.data(),
                                              static_cast<size_t>(g6_values.size()));
        return G6Calculator::fitCorrelationLength(r_span, g6_span);
      },
      py::arg("r_bins"), py::arg("g6_values"),
      "Fits the correlation length xi_6 from radial profile g6(r).\n\n"
      "Parameters\n----------\n"
      "r_bins : numpy.ndarray\n"
      "g6_values : numpy.ndarray\n\n"
      "Returns\n-------\n"
      "CorrelationLengthFit : Fitted correlation length xi_6, R^2, and plateau status.");
}
