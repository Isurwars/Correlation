#include "core_bindings.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace correlation::core;

void init_core(py::module_ &mod) {
  py::class_<Element>(mod, "Element")
      .def(py::init<>())
      .def_readwrite("symbol", &Element::symbol)
      .def_property(
          "id", [](const Element &element) { return element.id.value; },
          [](Element &element, int value) { element.id.value = value; });

  py::class_<Atom>(mod, "Atom")
      .def(py::init<>())
      .def_property(
          "position", [](const Atom &atom) { return atom.position().array(); },
          [](Atom &atom, const std::array<real_t, 3> &pos) {
            atom.setPosition(correlation::math::Vector3<real_t>(pos));
          })
      .def_property("id", &Atom::id, &Atom::setID)
      .def_property("element", &Atom::element, &Atom::setElement);

  py::class_<Cell>(mod, "Cell")
      .def(py::init<>())
      .def(
          "add_atom",
          [](Cell &cell, const std::string &symbol, const std::array<real_t, 3> &pos) -> Atom & {
            return cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos));
          },
          py::return_value_policy::reference)
      .def("get_volume", &Cell::volume)
      .def_property("energy", &Cell::getEnergy, &Cell::setEnergy)
      .def_property_readonly(
          "atoms", [](const Cell &cell) -> const std::vector<Atom> & { return cell.atoms(); },
          py::return_value_policy::reference_internal)
      .def(
          "get_positions",
          [](const Cell &cell) -> py::array_t<real_t> {
            py::module_ warnings = py::module_::import("warnings");
            warnings.attr("warn")("get_positions() is deprecated, use the zero-copy .positions property instead.",
                                  warnings.attr("DeprecationWarning"));
            const auto &atoms = cell.atoms();
            const size_t num_atoms = atoms.size();
            py::array_t<real_t> arr({static_cast<py::ssize_t>(num_atoms), py::ssize_t(3)});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < num_atoms; ++i) {
              const auto &position = atoms[i].position();
              buf(i, 0) = position[0];
              buf(i, 1) = position[1];
              buf(i, 2) = position[2];
            }
            return arr;
          },
          "Deprecated: Return all atom positions as a NumPy array. Use .positions instead.")
      .def_property_readonly(
          "positions",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &cell = obj.cast<const Cell &>();
            const auto &atoms = cell.atoms();
            if (atoms.empty()) {
              return py::array_t<real_t>();
            }

            ssize_t stride_row = sizeof(Atom);
            ssize_t stride_col = sizeof(real_t);
            ssize_t rows = static_cast<ssize_t>(atoms.size());
            ssize_t cols = 3;

            const real_t *ptr = atoms[0].position().begin();

            return py::array_t<real_t>({rows, cols}, {stride_row, stride_col}, ptr, obj);
          },
          "Zero-copy access to atom positions as a (N, 3) NumPy array.")
      .def_property_readonly(
          "velocities",
          [](py::object &obj) -> py::array_t<real_t> {
            const auto &cell = obj.cast<const Cell &>();
            const auto &atoms = cell.atoms();
            if (atoms.empty()) {
              return py::array_t<real_t>();
            }

            ssize_t stride_row = sizeof(Atom);
            ssize_t stride_col = sizeof(real_t);
            ssize_t rows = static_cast<ssize_t>(atoms.size());
            ssize_t cols = 3;

            const real_t *ptr = atoms[0].velocity().begin();

            return py::array_t<real_t>({rows, cols}, {stride_row, stride_col}, ptr, obj);
          },
          "Zero-copy access to atom velocities as a (N, 3) NumPy array.")
      .def(
          "get_element_ids",
          [](const Cell &cell) -> py::array_t<int> {
            const auto &atoms = cell.atoms();
            const size_t num_atoms = atoms.size();
            py::array_t<int> arr(static_cast<py::ssize_t>(num_atoms));
            auto buf = arr.mutable_unchecked<1>();
            for (size_t i = 0; i < num_atoms; ++i) {
              buf(i) = atoms[i].element_id();
            }
            return arr;
          },
          "Return element type IDs for all atoms as a NumPy array of shape (N,).")
      .def(
          "get_lattice_parameters",
          [](const Cell &cell) -> py::array_t<real_t> {
            const auto &lattice_parameters = cell.lattice_parameters();
            py::array_t<real_t> arr(6);
            auto buf = arr.mutable_unchecked<1>();
            for (int i = 0; i < 6; ++i) {
              buf(i) = lattice_parameters.at(i);
            }
            return arr;
          },
          "Return lattice parameters [a, b, c, alpha, beta, gamma] as a NumPy array.");

  py::class_<Trajectory>(mod, "Trajectory")
      .def(py::init<>())
      .def_property("time_step", &Trajectory::getTimeStep, &Trajectory::setTimeStep)
      .def("num_frames", &Trajectory::getFrameCount)
      .def("__len__", &Trajectory::getFrameCount)
      .def("__getitem__",
           [](const Trajectory &trajectory, int index) -> Cell {
             int count = static_cast<int>(trajectory.getFrameCount());
             if (index < 0) {
               index += count;
             }
             if (index < 0 || index >= count) {
               throw py::index_error("Trajectory index out of range");
             }
             return trajectory.getFrame(static_cast<size_t>(index));
           })
      .def_property_readonly(
          "frames", [](Trajectory &trajectory) -> std::vector<Cell> & { return trajectory.getFrames(); },
          py::return_value_policy::reference_internal);
}
