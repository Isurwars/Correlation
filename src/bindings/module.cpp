/**
 * @file module.cpp
 * @brief Root pybind11 module definition for the `correlation` Python package.
 *
 * Registers all sub-layers in dependency order:
 *   1. math        — KernelType enum
 *   2. core        — Atom, Cell, Trajectory
 *   3. io          — file reading (ReaderFactory)
 *   4. analysis    — AnalysisSettings, Histogram, StructureAnalyzer,
 *                    TrajectoryAnalyzer, DistributionFunctions
 *   5. calculators — BaseCalculator, CalculatorFactory
 *   6. writers     — BaseWriter, CSVWriter, WriterFactory, write_csv()
 */

#include "core_bindings.hpp"
#include "io_bindings.hpp"
#include "math_bindings.hpp"
#include "analysis_bindings.hpp"
#include "calculator_bindings.hpp"
#include "writer_bindings.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(correlation, m) {
    m.doc() =
        "Correlation: Liquid and Amorphous Solid Analysis Tool — Python Bindings\n\n"
        "Quick-start example::\n\n"
        "    import correlation\n\n"
        "    traj = correlation.read('my_sim.lammps')\n"
        "    cell = traj.frames[0]\n\n"
        "    df = correlation.DistributionFunctions(cell, cutoff=6.0)\n"
        "    df.calculate_rdf(r_max=15.0, bin_width=0.02)\n"
        "    df.smooth_all(sigma=0.05)\n\n"
        "    hist = df.get_histogram('g(r)')\n"
        "    print(hist.bins[:5], hist.partials['Total'][:5])\n\n"
        "    correlation.write_csv('output/sample', df)\n";

    init_math(m);
    init_core(m);
    init_io(m);
    init_analysis(m);
    init_calculators(m);
    init_writers(m);
}
