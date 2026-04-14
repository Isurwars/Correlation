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
#include "calculators/BaseCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>

namespace py = pybind11;
using namespace correlation::calculators;
using namespace correlation::analysis;

// ------------------------------------------------------------------
// Trampoline — lets users subclass BaseCalculator in Python
// ------------------------------------------------------------------
class PyBaseCalculator : public BaseCalculator {
public:
    using BaseCalculator::BaseCalculator;

    std::string getName() const override {
        PYBIND11_OVERRIDE_PURE(std::string, BaseCalculator, getName);
    }
    std::string getShortName() const override {
        PYBIND11_OVERRIDE_PURE(std::string, BaseCalculator, getShortName);
    }
    std::string getGroup() const override {
        PYBIND11_OVERRIDE_PURE(std::string, BaseCalculator, getGroup);
    }
    std::string getDescription() const override {
        PYBIND11_OVERRIDE_PURE(std::string, BaseCalculator, getDescription);
    }
    bool isFrameCalculator() const override {
        PYBIND11_OVERRIDE_PURE(bool, BaseCalculator, isFrameCalculator);
    }
    bool isTrajectoryCalculator() const override {
        PYBIND11_OVERRIDE_PURE(bool, BaseCalculator, isTrajectoryCalculator);
    }
    void calculateFrame(DistributionFunctions &df,
                        const AnalysisSettings &settings) const override {
        PYBIND11_OVERRIDE(void, BaseCalculator, calculateFrame, df, settings);
    }
    void calculateTrajectory(DistributionFunctions &df,
                             const correlation::core::Trajectory &traj,
                             const AnalysisSettings &settings) const override {
        PYBIND11_OVERRIDE(void, BaseCalculator, calculateTrajectory, df, traj, settings);
    }
};

void init_calculators(py::module_ &m) {
    // ------------------------------------------------------------------
    // BaseCalculator
    // ------------------------------------------------------------------
    py::class_<BaseCalculator, PyBaseCalculator>(m, "BaseCalculator",
        "Abstract base class for all analysis calculators.\n\n"
        "Subclass this in Python to create custom calculators.")
        .def(py::init<>())
        .def("get_name",        &BaseCalculator::getName,
             "Full display name of the calculator (e.g. 'g(r), J(r), G(r)').")
        .def("get_short_name",  &BaseCalculator::getShortName,
             "Short identifier used as a key (e.g. 'RDF').")
        .def("get_group",       &BaseCalculator::getGroup,
             "UI group this calculator belongs to (e.g. 'Radial', 'Angular').")
        .def("get_description", &BaseCalculator::getDescription,
             "Human-readable description of the calculator.")
        .def("is_frame_calculator",      &BaseCalculator::isFrameCalculator,
             "True if this calculator operates per-frame.")
        .def("is_trajectory_calculator", &BaseCalculator::isTrajectoryCalculator,
             "True if this calculator operates over the full trajectory.");

    // ------------------------------------------------------------------
    // CalculatorFactory access — module-level free functions
    // (CalculatorFactory holds unique_ptrs and cannot be copy-constructed)
    // ------------------------------------------------------------------
    m.def("list_calculators",
        []() {
            std::vector<std::string> names;
            for (const auto &c : CalculatorFactory::instance().getCalculators())
                names.push_back(c->getShortName());
            return names;
        },
        "Return a list of short names for all registered calculators\n"
        "(e.g. ['RDF', 'SQ', 'PAD', ...]).");

    m.def("get_calculator",
        [](const std::string &name) -> const BaseCalculator * {
            const auto *c = CalculatorFactory::instance().getCalculator(name);
            if (!c)
                throw std::runtime_error("No calculator registered with name: " + name);
            return c;
        },
        py::arg("name"),
        py::return_value_policy::reference,
        "Look up a registered calculator by its short name (e.g. 'RDF').\n\n"
        "Returns a reference into the factory's singleton registry.\n\n"
        "Parameters\n----------\n"
        "name : str\n"
        "    Short name of the desired calculator.\n\n"
        "Raises\n------\n"
        "RuntimeError\n"
        "    If no calculator with that name is registered.");

    m.def("get_all_calculators",
        []() {
            std::vector<const BaseCalculator *> out;
            for (const auto &c : CalculatorFactory::instance().getCalculators())
                out.push_back(c.get());
            return out;
        },
        py::return_value_policy::reference,
        "Return a list of all registered calculator objects.");
}
