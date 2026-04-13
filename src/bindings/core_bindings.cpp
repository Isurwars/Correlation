#include "core_bindings.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace correlation::core;

void init_core(py::module_ &m) {
    py::class_<Element>(m, "Element")
        .def(py::init<>())
        .def_readwrite("symbol", &Element::symbol)
        .def_property("id", [](const Element &e){ return e.id.value; }, [](Element &e, int v){ e.id.value = v; });

    py::class_<Atom>(m, "Atom")
        .def(py::init<>())
        .def_property("position", 
            [](const Atom &a) { return a.position().array(); },
            [](Atom &a, const std::array<double, 3> &pos) { a.setPosition(correlation::math::Vector3<double>(pos)); })
        .def_property("id", &Atom::id, &Atom::setID)
        .def_property("element", &Atom::element, &Atom::setElement);

    py::class_<Cell>(m, "Cell")
        .def(py::init<>())
        .def("add_atom", [](Cell &c, const std::string &symbol, const std::array<double, 3> &pos) -> Atom& {
            return c.addAtom(symbol, correlation::math::Vector3<double>(pos));
        }, py::return_value_policy::reference)
        .def("get_volume", &Cell::volume)
        .def_property("energy", &Cell::getEnergy, &Cell::setEnergy)
        .def_property_readonly("atoms", &Cell::atoms);

    py::class_<Trajectory>(m, "Trajectory")
        .def(py::init<>())
        .def_property("time_step", &Trajectory::getTimeStep, &Trajectory::setTimeStep)
        .def("num_frames", &Trajectory::getFrameCount)
        .def_property_readonly("frames", [](Trajectory &t) -> std::vector<Cell>& { return t.getFrames(); }, py::return_value_policy::reference_internal);
}
