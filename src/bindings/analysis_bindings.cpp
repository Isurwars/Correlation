/**
 * @file analysis_bindings.cpp
 * @brief Pybind11 bindings for correlation::analysis types.
 *
 * Exposes (in dependency order):
 *   - AnalysisSettings
 *   - Histogram
 *   - StructureAnalyzer
 *   - TrajectoryAnalyzer
 *   - DistributionFunctions  (including static compute_mean)
 */

#include "analysis_bindings.hpp"

#include "analysis/DistributionFunctions.hpp"
#include "analysis/StructureAnalyzer.hpp"
#include "analysis/TrajectoryAnalyzer.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Smoothing.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>

namespace py = pybind11;
using namespace correlation::analysis;
using namespace correlation::core;
using namespace correlation::math;

void init_analysis(py::module_ &m) {

    // ------------------------------------------------------------------
    // AnalysisSettings
    // ------------------------------------------------------------------
    py::class_<AnalysisSettings>(m, "AnalysisSettings",
        "Configuration bag for all distribution-function calculations.\n\n"
        "All parameters have sensible defaults so only the fields you want\n"
        "to override need to be set.")
        .def(py::init<>())
        // Radial
        .def_readwrite("r_max",            &AnalysisSettings::r_max,
                       "Maximum radius for RDF calculations (Å). Default 20.0.")
        .def_readwrite("r_bin_width",      &AnalysisSettings::r_bin_width,
                       "Bin width for radial distributions (Å). Default 0.02.")
        .def_readwrite("r_int_max",        &AnalysisSettings::r_int_max,
                       "Cutoff for integration-based properties (Å). Default 10.0.")
        // Reciprocal space
        .def_readwrite("q_max",            &AnalysisSettings::q_max,
                       "Maximum momentum transfer for S(Q) (Å⁻¹). Default 20.0.")
        .def_readwrite("q_bin_width",      &AnalysisSettings::q_bin_width,
                       "Bin width for S(Q) (Å⁻¹). Default 0.02.")
        // Angular
        .def_readwrite("angle_bin_width",    &AnalysisSettings::angle_bin_width,
                       "Bin width for bond angle distributions (°). Default 1.0.")
        .def_readwrite("dihedral_bin_width", &AnalysisSettings::dihedral_bin_width,
                       "Bin width for dihedral distributions (°). Default 1.0.")
        // Rings
        .def_readwrite("max_ring_size",    &AnalysisSettings::max_ring_size,
                       "Maximum ring size to search for. Default 8.")
        // Calculator selection
        .def_readwrite("active_calculators", &AnalysisSettings::active_calculators,
                       "Dict mapping calculator IDs to enabled state.\n"
                       "Empty dict means all calculators are enabled.")
        // Smoothing
        .def_readwrite("smoothing",        &AnalysisSettings::smoothing,
                       "Whether to apply post-processing smoothing. Default True.")
        .def_readwrite("smoothing_sigma",  &AnalysisSettings::smoothing_sigma,
                       "Gaussian smoothing standard deviation. Default 0.1.")
        .def_readwrite("smoothing_kernel", &AnalysisSettings::smoothing_kernel,
                       "Kernel type for smoothing (KernelType enum). Default Gaussian.")
        .def("is_active",
            [](const AnalysisSettings &s, const std::string &id) {
                return s.isActive(id);
            },
            py::arg("id"),
            "Return True if the given calculator ID is enabled.");

    // ------------------------------------------------------------------
    // Histogram
    // ------------------------------------------------------------------
    py::class_<Histogram>(m, "Histogram",
        "Container for a single calculated distribution function.\n\n"
        "Holds the x-axis bins, all partial distributions, and optional\n"
        "smoothed variants.")
        .def(py::init<>())
        .def_readwrite("bins",             &Histogram::bins,
                       "X-axis values (radii, angles, q-values, etc.).")
        .def_readwrite("title",            &Histogram::title,
                       "Descriptive plot title.")
        .def_readwrite("x_label",          &Histogram::x_label,
                       "X-axis label string.")
        .def_readwrite("y_label",          &Histogram::y_label,
                       "Y-axis label string.")
        .def_readwrite("x_unit",           &Histogram::x_unit,
                       "Physical unit for the x-axis (e.g. 'Å').")
        .def_readwrite("y_unit",           &Histogram::y_unit,
                       "Physical unit for the y-axis.")
        .def_readwrite("description",      &Histogram::description,
                       "Internal description of what this histogram represents.")
        .def_readwrite("file_suffix",      &Histogram::file_suffix,
                       "Default suffix used when saving this histogram to disk.")
        .def_readwrite("partials",         &Histogram::partials,
                       "Dict mapping partial key (e.g. 'Si-O') to y-values.")
        .def_readwrite("smoothed_partials",&Histogram::smoothed_partials,
                       "Dict mapping partial key to smoothed y-values.");

    // ------------------------------------------------------------------
    // StructureAnalyzer
    // ------------------------------------------------------------------
    py::class_<StructureAnalyzer>(m, "StructureAnalyzer",
        "Computes pairwise distances, bond angles, and dihedral angles\n"
        "for a single simulation cell (frame).\n\n"
        "The tensors are indexed by element type: distances[e1][e2][pair_idx],\n"
        "angles[center][e1][e2][angle_idx], dihedrals[e1][e2][e3][e4][idx].")
        .def(py::init<Cell &, double,
                      const std::vector<std::vector<double>> &, bool>(),
             py::arg("cell"),
             py::arg("cutoff"),
             py::arg("bond_cutoffs_sq"),
             py::arg("ignore_periodic_self_interactions") = true,
             "Construct and immediately compute all pair data.\n\n"
             "Parameters\n----------\n"
             "cell : Cell\n"
             "    The periodic simulation cell.\n"
             "cutoff : float\n"
             "    Neighbor search cutoff radius (Å).\n"
             "bond_cutoffs_sq : list[list[float]]\n"
             "    Per-element-pair squared bond cutoffs.\n"
             "ignore_periodic_self_interactions : bool\n"
             "    If True, atoms do not interact with their own periodic images.")
        .def("distances", &StructureAnalyzer::distances,
             py::return_value_policy::reference_internal,
             "3D tensor of pairwise distances [e1][e2][pair_idx].")
        .def("angles",    &StructureAnalyzer::angles,
             py::return_value_policy::reference_internal,
             "4D tensor of bond angles [center_e][e1][e2][angle_idx].")
        .def("dihedrals", &StructureAnalyzer::dihedrals,
             py::return_value_policy::reference_internal,
             "5D tensor of dihedral angles [e1][e2][e3][e4][dihedral_idx].");

    // ------------------------------------------------------------------
    // TrajectoryAnalyzer
    // ------------------------------------------------------------------
    py::class_<TrajectoryAnalyzer>(m, "TrajectoryAnalyzer",
        "Orchestrates structural analysis across multiple frames of a trajectory.\n\n"
        "Provides per-frame StructureAnalyzer factories and trajectory metadata.")
        .def(py::init<Trajectory &, double,
                      const std::vector<std::vector<double>> &,
                      size_t, long long, bool,
                      std::function<void(float, const std::string &)>>(),
             py::arg("trajectory"),
             py::arg("neighbor_cutoff"),
             py::arg("bond_cutoffs"),
             py::arg("start_frame") = 0,
             py::arg("end_frame")   = -1LL,
             py::arg("ignore_periodic_self_interactions") = true,
             py::arg("progress_callback") = nullptr,
             "Construct a TrajectoryAnalyzer.\n\n"
             "Parameters\n----------\n"
             "trajectory : Trajectory\n"
             "    The trajectory to analyze.\n"
             "neighbor_cutoff : float\n"
             "    Global neighbor search cutoff radius (Å).\n"
             "bond_cutoffs : list[list[float]]\n"
             "    Per-element-pair bond cutoffs (Å).\n"
             "start_frame : int, optional\n"
             "    Index of the first frame to analyze. Default 0.\n"
             "end_frame : int, optional\n"
             "    Index of the last frame (-1 for all). Default -1.\n"
             "ignore_periodic_self_interactions : bool, optional\n"
             "    Skip atom–own-image interactions. Default True.\n"
             "progress_callback : callable, optional\n"
             "    Called with (fraction: float, message: str) during computation.")
        .def("get_num_frames",  &TrajectoryAnalyzer::getNumFrames,
             "Total number of frames in the analysis window.")
        .def("get_start_frame", &TrajectoryAnalyzer::getStartFrame,
             "Index of the first frame being analyzed.")
        .def("get_time_step",   &TrajectoryAnalyzer::getTimeStep,
             "Time step between frames (from trajectory metadata).")
        .def("get_neighbor_cutoff", &TrajectoryAnalyzer::getNeighborCutoff,
             "Global neighbor search cutoff radius (Å).")
        .def("create_analyzer",
            [](const TrajectoryAnalyzer &ta, size_t frame_idx) {
                return ta.createAnalyzer(frame_idx);
            },
            py::arg("frame_idx"),
            "Create a StructureAnalyzer for the given frame index.");

    // ------------------------------------------------------------------
    // DistributionFunctions
    // ------------------------------------------------------------------
    py::class_<DistributionFunctions>(m, "DistributionFunctions",
        "Manager for distribution function calculations (RDF, PAD, S(Q), …).\n\n"
        "This is the primary analysis object. Construct it with a Cell (or obtain\n"
        "a trajectory-averaged instance via compute_mean), then call the desired\n"
        "calculate_* methods, and retrieve results via get_histogram().\n\n"
        ".. note::\n"
        "   The DistributionFunctions holds an internal reference to the Cell\n"
        "   passed at construction. The Cell (and owning Trajectory) must remain\n"
        "   alive for the lifetime of this object.")
        .def(py::init<Cell &, double,
                      const std::vector<std::vector<double>> &>(),
             py::arg("cell"),
             py::arg("cutoff") = 0.0,
             py::arg("bond_cutoffs") = std::vector<std::vector<double>>{},
             "Construct a DistributionFunctions for a single Cell.\n\n"
             "Parameters\n----------\n"
             "cell : Cell\n"
             "    The simulation cell to analyze.\n"
             "cutoff : float, optional\n"
             "    Neighbor search cutoff (Å). 0 uses a heuristic. Default 0.0.\n"
             "bond_cutoffs : list[list[float]], optional\n"
             "    Per-element-pair bond cutoffs. Default is empty (auto).")

        // Accessors
        .def("get_histogram",
            [](const DistributionFunctions &df, const std::string &name) -> const Histogram & {
                return df.getHistogram(name);
            },
            py::arg("name"),
            py::return_value_policy::reference_internal,
            "Return the Histogram for the given key (e.g. 'g(r)', 'S(Q)').\n"
            "Raises KeyError if the histogram has not been calculated.")
        .def("get_all_histograms",
            [](const DistributionFunctions &df) -> const std::map<std::string, Histogram> & {
                return df.getAllHistograms();
            },
            py::return_value_policy::reference_internal,
            "Return a dict of all calculated histograms.")
        .def("get_available_histograms", &DistributionFunctions::getAvailableHistograms,
             "Return a list of names for all currently available histograms.")
        .def("get_ashcroft_weights",
            [](const DistributionFunctions &df) -> const std::map<std::string, double> & {
                return df.getAshcroftWeights();
            },
            py::return_value_policy::reference_internal,
            "Return the Ashcroft-Langreth weights used for S(Q) partials.")

        // Calculations
        .def("calculate_rdf",
            &DistributionFunctions::calculateRDF,
            py::arg("r_max") = 20.0, py::arg("bin_width") = 0.05,
            "Calculate g(r), J(r), and G(r).\n\n"
            "Parameters\n----------\n"
            "r_max : float\n    Maximum radius (Å). Default 20.0.\n"
            "bin_width : float\n    Histogram bin width (Å). Default 0.05.")
        .def("calculate_pad",
            &DistributionFunctions::calculatePAD,
            py::arg("bin_width") = 1.0,
            "Calculate the Plane Angle Distribution (PAD).\n\n"
            "bin_width : float\n    Angular bin width (°). Default 1.0.")
        .def("calculate_dad",
            &DistributionFunctions::calculateDAD,
            py::arg("bin_width") = 1.0,
            "Calculate the Dihedral Angle Distribution (DAD).\n\n"
            "bin_width : float\n    Angular bin width (°). Default 1.0.")
        .def("calculate_vacf",
            [](DistributionFunctions &df, const Trajectory &traj,
               int max_correlation_frames, size_t start_frame, size_t end_frame) {
                df.calculateVACF(traj, max_correlation_frames,
                                 start_frame, end_frame);
            },
            py::arg("traj"),
            py::arg("max_correlation_frames") = -1,
            py::arg("start_frame") = 0,
            py::arg("end_frame")   = static_cast<size_t>(-1),
            "Calculate the Velocity Autocorrelation Function (VACF).\n\n"
            "Parameters\n----------\n"
            "traj : Trajectory\n    Trajectory with pre-calculated velocities.\n"
            "max_correlation_frames : int\n    Max lag frames (-1 = half trajectory).\n"
            "start_frame : int\n    First frame to include. Default 0.\n"
            "end_frame : int\n    Last frame (exclusive). Default all.")
        .def("calculate_vdos",
            &DistributionFunctions::calculateVDOS,
            "Calculate the Vibrational Density of States (VDOS) from the VACF.\n"
            "Requires calculate_vacf() to have been called first.")
        .def("calculate_xrd",
            &DistributionFunctions::calculateXRD,
            py::arg("lambda")    = 1.5406,
            py::arg("theta_min") = 5.0,
            py::arg("theta_max") = 90.0,
            py::arg("bin_width") = 1.0,
            "Calculate the X-Ray Diffraction (XRD) pattern.\n\n"
            "Parameters\n----------\n"
            "lambda : float\n    X-ray wavelength (Å). Default 1.5406 (Cu Kα).\n"
            "theta_min : float\n    Minimum 2θ angle (°). Default 5.0.\n"
            "theta_max : float\n    Maximum 2θ angle (°). Default 90.0.\n"
            "bin_width : float\n    2θ bin width (°). Default 1.0.")
        .def("calculate_cn",
            &DistributionFunctions::calculateCoordinationNumber,
            "Calculate the Coordination Number (CN) distribution.\n"
            "Requires neighbors to have been computed (non-zero cutoff).")

        // Smoothing
        .def("smooth",
            [](DistributionFunctions &df, const std::string &name, double sigma,
               KernelType kernel) {
                df.smooth(name, sigma, kernel);
            },
            py::arg("name"), py::arg("sigma"),
            py::arg("kernel") = KernelType::Gaussian,
            "Smooth a specific histogram.\n\n"
            "Parameters\n----------\n"
            "name : str\n    Histogram name (e.g. 'g(r)').\n"
            "sigma : float\n    Kernel bandwidth.\n"
            "kernel : KernelType\n    Smoothing kernel. Default Gaussian.")
        .def("smooth_all",
            [](DistributionFunctions &df, double sigma, KernelType kernel) {
                df.smoothAll(sigma, kernel);
            },
            py::arg("sigma"),
            py::arg("kernel") = KernelType::Gaussian,
            "Smooth all available histograms.\n\n"
            "sigma : float\n    Kernel bandwidth.\n"
            "kernel : KernelType\n    Smoothing kernel. Default Gaussian.")

        // Accumulation helpers
        .def("add",   &DistributionFunctions::add,   py::arg("other"),
             "Accumulate histograms from another DistributionFunctions object\n"
             "(used for trajectory averaging).")
        .def("scale", &DistributionFunctions::scale, py::arg("factor"),
             "Scale all histogram values by *factor* (used for normalization).")

        // Static trajectory-mean method
        .def_static("compute_mean",
            [](Trajectory &trajectory,
               const TrajectoryAnalyzer &analyzer,
               size_t start_frame,
               const AnalysisSettings &settings,
               std::function<void(float, const std::string &)> progress_callback) {
                return DistributionFunctions::computeMean(
                    trajectory, analyzer, start_frame, settings,
                    progress_callback);
            },
            py::arg("trajectory"),
            py::arg("analyzer"),
            py::arg("start_frame") = 0,
            py::arg("settings")    = AnalysisSettings{},
            py::arg("progress_callback") = nullptr,
            "Compute trajectory-averaged distribution functions in parallel.\n\n"
            "Parameters\n----------\n"
            "trajectory : Trajectory\n"
            "    The trajectory to analyze.\n"
            "analyzer : TrajectoryAnalyzer\n"
            "    Pre-configured analyzer providing frame-level neighbor info.\n"
            "start_frame : int, optional\n"
            "    Frame to start from. Default 0.\n"
            "settings : AnalysisSettings, optional\n"
            "    Analysis configuration (bin widths, cutoffs, enabled calcs).\n"
            "progress_callback : callable, optional\n"
            "    Called with (fraction: float, message: str) during computation.\n\n"
            "Returns\n-------\n"
            "DistributionFunctions\n"
            "    A new object containing the averaged results.");
}
