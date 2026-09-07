#include "core_bindings.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using correlation::real_t;
using namespace correlation::core;

void init_core(py::module_ &mod) {
  py::class_<Element>(mod, "Element", "Chemical element specification.")
      .def(py::init<>(), "Construct an empty Element.")
      .def_readwrite("symbol", &Element::symbol, "Chemical symbol (e.g., 'C', 'O', 'Fe').")
      .def_property(
          "id", [](const Element &element) { return element.id.value; },
          [](Element &element, int value) { element.id.value = value; },
          "Unique integer identifier for the element type.");

  py::class_<Atom>(mod, "Atom", "Single atomic site with Cartesian position and chemical element.")
      .def(py::init<>(), "Construct an atom at the origin.")
      .def_property(
          "position", [](const Atom &atom) { return atom.position().array(); },
          [](Atom &atom, const std::array<real_t, 3> &pos) {
            atom.setPosition(correlation::math::Vector3<real_t>(pos));
          },
          "Cartesian position [x, y, z] in Angstroms.")
      .def_property("id", &Atom::id, &Atom::setID, "Atom identifier.")
      .def_property("element", &Atom::element, &Atom::setElement, "Chemical element metadata.");

  py::class_<Cell>(mod, "Cell", "Simulation cell containing atomic coordinates and periodic box geometry.")
      .def(py::init<>(), "Construct an empty non-periodic Cell.")
      .def(py::init<const std::array<real_t, 6> &>(), py::arg("lattice_parameters"),
           "Construct a Cell from lattice parameters [a, b, c, alpha, beta, gamma].")
      .def(py::init([](const std::array<real_t, 3> &vec_a, const std::array<real_t, 3> &vec_b,
                       const std::array<real_t, 3> &vec_c) {
             return Cell(correlation::math::Vector3<real_t>(vec_a), correlation::math::Vector3<real_t>(vec_b),
                         correlation::math::Vector3<real_t>(vec_c));
           }),
           py::arg("a"), py::arg("b"), py::arg("c"), "Construct a Cell from three lattice vectors a, b, and c.")
      .def(
          "add_atom",
          [](Cell &cell, const std::string &symbol, const std::array<real_t, 3> &pos) -> Atom & {
            return cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos));
          },
          py::arg("symbol"), py::arg("position"), py::return_value_policy::reference,
          "Add an atom with given symbol and Cartesian position to the cell.")
      .def("get_volume", &Cell::volume, "Get unit cell volume in cubic Angstroms.")
      .def_property_readonly("volume", &Cell::volume, "Unit cell volume in cubic Angstroms.")
      .def_property("energy", &Cell::getEnergy, &Cell::setEnergy, "Potential energy of the cell snapshot.")
      .def_property_readonly("atom_count", &Cell::atomCount, "Number of atoms in the cell.")
      .def("__len__", &Cell::atomCount, "Number of atoms in the cell.")
      .def(
          "__iter__", [](const Cell &cell) { return py::make_iterator(cell.atoms().begin(), cell.atoms().end()); },
          py::keep_alive<0, 1>(), "Iterate over atoms in the cell.")
      .def_property_readonly(
          "atoms", [](const Cell &cell) -> const std::vector<Atom> & { return cell.atoms(); },
          py::return_value_policy::reference_internal, "List of atoms contained in this cell.")
      .def(
          "get_positions",
          [](const Cell &cell) -> py::array_t<real_t> {
            const py::module_ warnings = py::module_::import("warnings");
            warnings.attr("warn")("get_positions() is deprecated, use the zero-copy .positions property instead.",
                                  warnings.attr("DeprecationWarning"));
            const auto &atoms = cell.atoms();
            const size_t num_atoms = atoms.size();
            py::array_t<real_t> arr({static_cast<py::ssize_t>(num_atoms), static_cast<py::ssize_t>(3)});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t idx = 0; idx < num_atoms; ++idx) {
              const auto &position = atoms[idx].position();
              buf(idx, 0) = position[0];
              buf(idx, 1) = position[1];
              buf(idx, 2) = position[2];
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
              return {};
            }

            const py::ssize_t stride_row = sizeof(Atom);
            const py::ssize_t stride_col = sizeof(real_t);
            const auto rows = static_cast<py::ssize_t>(atoms.size());
            const py::ssize_t cols = 3;

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
              return {};
            }

            const py::ssize_t stride_row = sizeof(Atom);
            const py::ssize_t stride_col = sizeof(real_t);
            const auto rows = static_cast<py::ssize_t>(atoms.size());
            const py::ssize_t cols = 3;

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
            for (size_t idx = 0; idx < num_atoms; ++idx) {
              buf(idx) = atoms[idx].element_id();
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
            for (int idx = 0; idx < 6; ++idx) {
              buf(idx) = lattice_parameters.at(idx);
            }
            return arr;
          },
          "Return lattice parameters [a, b, c, alpha, beta, gamma] as a NumPy array.")
      .def_property_readonly(
          "lattice_vectors",
          [](const Cell &cell) -> py::array_t<real_t> {
            const auto &lattice_matrix = cell.latticeVectors();
            py::array_t<real_t> arr({static_cast<py::ssize_t>(3), static_cast<py::ssize_t>(3)});
            auto buf = arr.mutable_unchecked<2>();
            for (py::ssize_t row_idx = 0; row_idx < 3; ++row_idx) {
              for (py::ssize_t col_idx = 0; col_idx < 3; ++col_idx) {
                buf(row_idx, col_idx) = lattice_matrix[row_idx][col_idx];
              }
            }
            return arr;
          },
          "Lattice vectors as a (3, 3) NumPy array where rows are vectors a, b, and c.");

  py::class_<Trajectory>(mod, "Trajectory", "Time-series collection of Cell simulation snapshots.")
      .def(py::init<>(), "Construct an empty trajectory.")
      .def_property("time_step", &Trajectory::getTimeStep, &Trajectory::setTimeStep,
                    "Time interval between consecutive frames.")
      .def("num_frames", &Trajectory::getFrameCount, "Number of frames in the trajectory.")
      .def("__len__", &Trajectory::getFrameCount, "Number of frames in the trajectory.")
      .def(
          "__iter__",
          [](Trajectory &trajectory) {
            auto &frames = trajectory.getFrames();
            return py::make_iterator(frames.begin(), frames.end());
          },
          py::keep_alive<0, 1>(), "Iterate over frames in the trajectory.")
      .def(
          "__getitem__",
          [](const Trajectory &trajectory, int index) -> Cell {
            const int count = static_cast<int>(trajectory.getFrameCount());
            if (index < 0) {
              index += count;
            }
            if (index < 0 || index >= count) {
              throw py::index_error("Trajectory index out of range");
            }
            return trajectory.getFrame(static_cast<size_t>(index));
          },
          py::arg("index"), "Access a frame snapshot by integer index.")
      .def("add_frame", &Trajectory::addFrame, py::arg("frame"), "Append a Cell frame to the trajectory.")
      .def("append", &Trajectory::addFrame, py::arg("frame"),
           "Append a Cell frame to the trajectory (alias for add_frame).")
      .def_property_readonly(
          "frames", [](Trajectory &trajectory) -> std::vector<Cell> & { return trajectory.getFrames(); },
          py::return_value_policy::reference_internal, "Direct reference to the sequence of Cell frames.");
}
