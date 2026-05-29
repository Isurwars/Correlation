#include "core_bindings.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace correlation::core;

void init_core(py::module_ &m) {
  py::class_<Element>(m, "Element")
      .def(py::init<>())
      .def_readwrite("symbol", &Element::symbol)
      .def_property("id", [](const Element &e) { return e.id.value; }, [](Element &e, int v) { e.id.value = v; });

  py::class_<Atom>(m, "Atom")
      .def(py::init<>())
      .def_property(
          "position", [](const Atom &a) { return a.position().array(); },
          [](Atom &a, const std::array<double, 3> &pos) { a.setPosition(correlation::math::Vector3<double>(pos)); })
      .def_property("id", &Atom::id, &Atom::setID)
      .def_property("element", &Atom::element, &Atom::setElement);

  py::class_<Cell>(m, "Cell")
      .def(py::init<>())
      .def(
          "add_atom",
          [](Cell &c, const std::string &symbol, const std::array<double, 3> &pos) -> Atom & {
            return c.addAtom(symbol, correlation::math::Vector3<double>(pos));
          },
          py::return_value_policy::reference)
      .def("get_volume", &Cell::volume)
      .def_property("energy", &Cell::getEnergy, &Cell::setEnergy)
      .def_property_readonly(
          "atoms", [](const Cell &c) -> const std::vector<Atom> & { return c.atoms(); },
          py::return_value_policy::reference_internal)
      .def(
          "get_positions",
          [](const Cell &c) -> py::array_t<double> {
            py::module_ warnings = py::module_::import("warnings");
            warnings.attr("warn")("get_positions() is deprecated, use the zero-copy .positions property instead.",
                                  warnings.attr("DeprecationWarning"));
            const auto &atoms = c.atoms();
            const size_t n = atoms.size();
            py::array_t<double> arr({n, size_t(3)});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n; ++i) {
              const auto &p = atoms[i].position();
              buf(i, 0) = p[0];
              buf(i, 1) = p[1];
              buf(i, 2) = p[2];
            }
            return arr;
          },
          "Deprecated: Return all atom positions as a NumPy array. Use .positions instead.")
      .def_property_readonly(
          "positions",
          [](py::object &obj) -> py::array_t<double> {
            auto &c = obj.cast<const Cell &>();
            const auto &atoms = c.atoms();
            if (atoms.empty())
              return py::array_t<double>();

            ssize_t stride_row = sizeof(Atom);
            ssize_t stride_col = sizeof(double);
            ssize_t rows = atoms.size();
            ssize_t cols = 3;

            const double *ptr = atoms[0].position().begin();

            return py::array_t<double>({rows, cols}, {stride_row, stride_col}, ptr, obj);
          },
          "Zero-copy access to atom positions as a (N, 3) NumPy array.")
      .def_property_readonly(
          "velocities",
          [](py::object &obj) -> py::array_t<double> {
            auto &c = obj.cast<const Cell &>();
            const auto &atoms = c.atoms();
            if (atoms.empty())
              return py::array_t<double>();

            ssize_t stride_row = sizeof(Atom);
            ssize_t stride_col = sizeof(double);
            ssize_t rows = atoms.size();
            ssize_t cols = 3;

            const double *ptr = atoms[0].velocity().begin();

            return py::array_t<double>({rows, cols}, {stride_row, stride_col}, ptr, obj);
          },
          "Zero-copy access to atom velocities as a (N, 3) NumPy array.")
      .def(
          "get_element_ids",
          [](const Cell &c) -> py::array_t<int> {
            const auto &atoms = c.atoms();
            const size_t n = atoms.size();
            py::array_t<int> arr(n);
            auto buf = arr.mutable_unchecked<1>();
            for (size_t i = 0; i < n; ++i) {
              buf(i) = atoms[i].element_id();
            }
            return arr;
          },
          "Return element type IDs for all atoms as a NumPy array of shape (N,).")
      .def(
          "get_lattice_parameters",
          [](const Cell &c) -> py::array_t<double> {
            const auto &lp = c.lattice_parameters();
            py::array_t<double> arr(6);
            auto buf = arr.mutable_unchecked<1>();
            for (int i = 0; i < 6; ++i)
              buf(i) = lp[i];
            return arr;
          },
          "Return lattice parameters [a, b, c, alpha, beta, gamma] as a NumPy array.");

  py::class_<Trajectory>(m, "Trajectory")
      .def(py::init<>())
      .def_property("time_step", &Trajectory::getTimeStep, &Trajectory::setTimeStep)
      .def("num_frames", &Trajectory::getFrameCount)
      .def_property_readonly(
          "frames", [](Trajectory &t) -> std::vector<Cell> & { return t.getFrames(); },
          py::return_value_policy::reference_internal);
}
