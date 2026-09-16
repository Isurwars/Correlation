"""
Correlation: Liquid and Amorphous Solid Analysis Tool — Python Bindings

Quick-start example::

    import correlation

    traj = correlation.read('my_sim.lammps')
    cell = traj.frames[0]

    dists = correlation.DistributionFunctions(cell, cutoff=6.0)
    dists.calculate_rdf(r_max=15.0, bin_width=0.02)
    dists.smooth_all(sigma=0.05)

    hist = dists.get_histogram('g(r)')
    print(hist.bins[:5], hist.partials['Total'][:5])

    correlation.write_csv('output/sample', dists)
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
__all__: list[str] = ['AnalysisSettings', 'Atom', 'BCC', 'BaseCalculator', 'BaseWriter', 'Biweight', 'BondCutoffRange', 'Bump', 'CNALabel', 'CSVWriter', 'Cell', 'Cosine', 'DistributionFunctions', 'Element', 'Epanechnikov', 'FCC', 'Gaussian', 'GaussianRBFConfig', 'GraphDescriptors', 'HCP', 'Histogram', 'ICO', 'KernelType', 'MLIPCalculator', 'MLIPInterface', 'MLIPOutput', 'MotifProjectedTDOS', 'OTHER', 'OrbDescriptorConfig', 'PeriodicGraphBuilder', 'PeriodicGraphData', 'StructuralElectronicCorrelation', 'StructureAnalyzer', 'TDOSCalculator', 'TDOSParams', 'Trajectory', 'TrajectoryAnalyzer', 'Triweight', 'build_orb_graph', 'build_periodic_graph', 'compute_cna_descriptor', 'compute_coordination_embedding', 'compute_graph_spectrum', 'compute_ring_statistics_descriptor', 'correlate_cna', 'correlate_steinhardt', 'get_all_calculators', 'get_calculator', 'get_writer', 'get_writer_for_extension', 'list_calculators', 'list_writers', 'populate_descriptors', 'read', 'write_csv']
class AnalysisSettings:
    """
    Configuration bag for all distribution-function calculations.
    
    All parameters have sensible defaults so only the fields you want
    to override need to be set.
    """
    def __init__(self) -> None:
        ...
    def is_active(self, idx: str) -> bool:
        """
        Return True if the given calculator index is enabled.
        """
    @property
    def active_calculators(self) -> dict[str, bool]:
        """
        Dict mapping calculator IDs to enabled state.
        Empty dict means all calculators are enabled.
        """
    @active_calculators.setter
    def active_calculators(self, arg0: collections.abc.Mapping[str, bool]) -> None:
        ...
    @property
    def angle_bin_width(self) -> float:
        """
        Bin width for bond angle distributions (°). Default 1.0.
        """
    @angle_bin_width.setter
    def angle_bin_width(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def dihedral_bin_width(self) -> float:
        """
        Bin width for dihedral distributions (°). Default 1.0.
        """
    @dihedral_bin_width.setter
    def dihedral_bin_width(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def hyperuniformity_samples(self) -> int:
        """
        Number of random sampling points for hyperuniformity. Default 10000.
        """
    @hyperuniformity_samples.setter
    def hyperuniformity_samples(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def lef_cutoff(self) -> float:
        """
        Cutoff radius for local entropy integration (Å). Default 5.0.
        """
    @lef_cutoff.setter
    def lef_cutoff(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def lef_sigma(self) -> float:
        """
        Gaussian sigma for local entropy smoothing (Å). Default 0.2.
        """
    @lef_sigma.setter
    def lef_sigma(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def max_ring_size(self) -> int:
        """
        Maximum ring size to search for. Default 8.
        """
    @max_ring_size.setter
    def max_ring_size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def q_bin_width(self) -> float:
        """
        Bin width for S(Q) (Å⁻¹). Default 0.02.
        """
    @q_bin_width.setter
    def q_bin_width(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def q_max(self) -> float:
        """
        Maximum momentum transfer for S(Q) (Å⁻¹). Default 20.0.
        """
    @q_max.setter
    def q_max(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def r_bin_width(self) -> float:
        """
        Bin width for radial distributions (Å). Default 0.02.
        """
    @r_bin_width.setter
    def r_bin_width(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def r_int_max(self) -> float:
        """
        Cutoff for integration-based properties (Å). Default 10.0.
        """
    @r_int_max.setter
    def r_int_max(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def r_max(self) -> float:
        """
        Maximum radius for RDF calculations (Å). Default 20.0.
        """
    @r_max.setter
    def r_max(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def smoothing(self) -> bool:
        """
        Whether to apply post-processing smoothing. Default True.
        """
    @smoothing.setter
    def smoothing(self, arg0: bool) -> None:
        ...
    @property
    def smoothing_kernel(self) -> KernelType:
        """
        Kernel used for smoothing (KernelType enum). Default KernelType.Gaussian.
        """
    @smoothing_kernel.setter
    def smoothing_kernel(self, arg0: KernelType) -> None:
        ...
    @property
    def smoothing_sigma(self) -> float:
        """
        Gaussian smoothing standard deviation. Default 0.1.
        """
    @smoothing_sigma.setter
    def smoothing_sigma(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class Atom:
    """
    Single atomic site with Cartesian position and chemical element.
    """
    def __init__(self) -> None:
        """
        Construct an atom at the origin.
        """
    @property
    def element(self) -> Element:
        """
        Chemical element metadata.
        """
    @element.setter
    def element(self, arg1: Element) -> None:
        ...
    @property
    def id(self) -> int:
        """
        Atom identifier.
        """
    @id.setter
    def id(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def position(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        """
        Cartesian position [x, y, z] in Angstroms.
        """
    @position.setter
    def position(self, arg1: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> None:
        ...
class BaseCalculator:
    """
    Abstract base class for all analysis calculators.
    
    Subclass this in Python to create custom calculators.
    """
    def __init__(self) -> None:
        ...
    def get_description(self) -> str:
        """
        Human-readable description of the calculator.
        """
    def get_group(self) -> str:
        """
        UI group this calculator belongs to (e.g. 'Radial', 'Angular').
        """
    def get_name(self) -> str:
        """
        Full display name of the calculator (e.g. 'g(r), J(r), G(r)').
        """
    def get_short_name(self) -> str:
        """
        Short identifier used as a key (e.g. 'RDF').
        """
    def is_configured(self) -> bool:
        """
        True if the calculator has all required dependencies and models configured to execute.
        """
    def is_frame_calculator(self) -> bool:
        """
        True if this calculator operates per-frame.
        """
    def is_trajectory_calculator(self) -> bool:
        """
        True if this calculator operates over the full trajectory.
        """
class BaseWriter:
    """
    Abstract base class for all file-format writers.
    """
    def get_extensions(self) -> list[str]:
        """
        List of supported file extensions (e.g. ['.csv']).
        """
    def get_name(self) -> str:
        """
        Display name of this writer (e.g. 'CSV', 'HDF5').
        """
    def write(self, base_path: str, dists: DistributionFunctions, smoothing: bool = False) -> None:
        """
        Write distribution function data to file(s).
        
        Parameters
        ----------
        base_path : str
            Base name for output files (without extension).
        dists : DistributionFunctions
            The analysis results to write.
        smoothing : bool, optional
            If True, also write smoothed data. Default is False.
        """
class BondCutoffRange:
    """
    Squared minimum and maximum bond distance cutoff bounds.
    """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, min_sq: typing.SupportsFloat | typing.SupportsIndex, max_sq: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def max_sq(self) -> float:
        ...
    @max_sq.setter
    def max_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def min_sq(self) -> float:
        ...
    @min_sq.setter
    def min_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class CNALabel:
    """
    Canonical Common Neighbor Analysis structural motifs.
    
    Members:
    
      OTHER
    
      FCC
    
      HCP
    
      BCC
    
      ICO
    """
    BCC: typing.ClassVar[CNALabel]  # value = <CNALabel.BCC: 3>
    FCC: typing.ClassVar[CNALabel]  # value = <CNALabel.FCC: 1>
    HCP: typing.ClassVar[CNALabel]  # value = <CNALabel.HCP: 2>
    ICO: typing.ClassVar[CNALabel]  # value = <CNALabel.ICO: 4>
    OTHER: typing.ClassVar[CNALabel]  # value = <CNALabel.OTHER: 0>
    __members__: typing.ClassVar[dict[str, CNALabel]]  # value = {'OTHER': <CNALabel.OTHER: 0>, 'FCC': <CNALabel.FCC: 1>, 'HCP': <CNALabel.HCP: 2>, 'BCC': <CNALabel.BCC: 3>, 'ICO': <CNALabel.ICO: 4>}
    @typing.overload
    def __eq__(self, other: CNALabel) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: CNALabel) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class CSVWriter(BaseWriter):
    """
    Writes distribution function histograms to CSV files.
    
    For each histogram (e.g. g(r)) a separate .csv file is created
    with all partials as columns.
    """
    def __init__(self) -> None:
        ...
    def write_all_csvs(self, base_path: str, dists: DistributionFunctions, write_smoothed: bool = False) -> None:
        """
        Write all available histograms to individual CSV files.
        """
class Cell:
    """
    Simulation cell containing atomic coordinates and periodic box geometry.
    """
    @staticmethod
    def from_ase(atoms: typing.Any) -> typing.Any:
        """
        
        Convert an ASE Atoms object into a Correlation Cell.
        
        Parameters
        ----------
        atoms : ase.Atoms
            Source ASE atomic configuration.
        
        Returns
        -------
        correlation.Cell
            Constructed Simulation Cell with matching lattice vectors and atomic positions.
        
        Raises
        ------
        ImportError
            If ASE is not installed in the current environment.
        """
    @staticmethod
    def from_pymatgen(structure_or_molecule: typing.Any) -> typing.Any:
        """
        
        Convert a Pymatgen Structure (periodic) or Molecule (non-periodic) into a Correlation Cell.
        
        Parameters
        ----------
        structure_or_molecule : pymatgen.core.Structure or pymatgen.core.Molecule
            Source Pymatgen representation.
        
        Returns
        -------
        correlation.Cell
            Constructed Correlation Cell.
        
        Raises
        ------
        ImportError
            If Pymatgen is not installed.
        """
    def to_ase(self) -> typing.Any:
        """
        
        Convert a Correlation Cell into an ASE Atoms object.
        
        Returns
        -------
        ase.Atoms
            Converted ASE Atoms object with matching positions, symbols, and cell parameters.
        
        Raises
        ------
        ImportError
            If ASE is not installed in the current environment.
        """
    def to_atom_graphs(self, config: typing.Any = None, device: typing.Any = None, output_dtype: typing.Any = None, max_num_neighbors: int | None = None) -> typing.Any:
        """
        
        Convert a Correlation Cell or PeriodicGraphData into an orb_models AtomGraphs batch.
        
        Parameters
        ----------
        config : correlation.OrbDescriptorConfig, optional
            Configuration parameters for ORB graph extraction.
        device : torch.device or str, optional
            Target device for PyTorch tensors.
        output_dtype : torch.dtype, optional
            Floating-point dtype for tensors (e.g. torch.float32).
        max_num_neighbors : int, optional
            Cap on neighbors per node.
        
        Returns
        -------
        orb_models.common.atoms.batch.graph_batch.AtomGraphs
            Batched atom graphs on the target device.
        """
    def to_pymatgen(self) -> typing.Any:
        """
        
        Convert a Correlation Cell into a Pymatgen Structure (if periodic) or Molecule (if non-periodic).
        
        ----------
        cell : correlation.Cell
            Source Correlation simulation cell.
        
        Returns
        -------
        pymatgen.core.Structure or pymatgen.core.Molecule
            Converted Pymatgen data structure.
        
        Raises
        ------
        ImportError
            If Pymatgen is not installed.
        """
    def __init__(self, *args, **kwargs):
        ...
    def __iter__(self) -> collections.abc.Iterator[Atom]:
        """
        Iterate over atoms in the cell.
        """
    def __len__(self) -> int:
        """
        Number of atoms in the cell.
        """
    def add_atom(self, *args, **kwargs):
        ...
    def get_element_ids(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        Return element type IDs for all atoms as a NumPy array of shape (N,).
        """
    def get_lattice_parameters(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Return lattice parameters [a, b, c, alpha, beta, gamma] as a NumPy array.
        """
    def get_positions(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Deprecated: Return all atom positions as a NumPy array. Use .positions instead.
        """
    def get_volume(self) -> float:
        """
        Get unit cell volume in cubic Angstroms.
        """
    @property
    def atom_count(self):
        ...
    @property
    def atoms(self) -> list[Atom]:
        """
        List of atoms contained in this cell.
        """
    @property
    def energy(self) -> float:
        """
        Potential energy of the cell snapshot.
        """
    @energy.setter
    def energy(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def lattice_vectors(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Lattice vectors as a (3, 3) NumPy array where rows are vectors a, b, and c.
        """
    @property
    def positions(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to atom positions as a (N, 3) NumPy array.
        """
    @property
    def velocities(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to atom velocities as a (N, 3) NumPy array.
        """
    @property
    def volume(self) -> float:
        """
        Unit cell volume in cubic Angstroms.
        """
class DistributionFunctions:
    """
    Manager for distribution function calculations (RDF, PAD, S(Q), …).
    
    This is the primary analysis object. Construct it with a Cell (or obtain
    a trajectory-averaged instance via compute_mean), then call the desired
    calculate_* methods, and retrieve results via get_histogram().
    
    .. note::
       The DistributionFunctions holds an internal reference to the Cell
       passed at construction. The Cell (and owning Trajectory) must remain
       alive for the lifetime of this object.
    """
    @staticmethod
    def compute_mean(trajectory: Trajectory, analyzer: TrajectoryAnalyzer, start_frame: typing.SupportsInt | typing.SupportsIndex = 0, settings: AnalysisSettings = ..., progress_callback: collections.abc.Callable[[float, str], None] | None = None) -> DistributionFunctions:
        """
        Compute trajectory-averaged distribution functions in parallel.
        
        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to analyze.
        analyzer : TrajectoryAnalyzer
            Pre-configured analyzer providing frame-level neighbor info.
        start_frame : int, optional
            Frame to start from. Default 0.
        settings : AnalysisSettings, optional
            Analysis configuration (bin widths, cutoffs, enabled calcs).
        progress_callback : callable, optional
            Called with (fraction: float, message: str) during computation.
        
        Returns
        -------
        DistributionFunctions
            A new object containing the averaged results.
        """
    def __init__(self, cell, cutoff: float = 0.0, bond_cutoffs = None, **kwargs):
        ...
    def add(self, other: DistributionFunctions) -> None:
        """
        Accumulate histograms from another DistributionFunctions object
        (used for trajectory averaging).
        """
    def calculate_cn(self) -> None:
        """
        Calculate the Coordination Number (CN) distribution.
        Requires neighbors to have been computed (non-zero cutoff).
        """
    def calculate_dad(self, bin_width: typing.SupportsFloat | typing.SupportsIndex = 0.25) -> None:
        """
        Calculate the Dihedral Angle Distribution (DAD).
        
        bin_width : float
            Angular bin width (°). Default 0.25.
        """
    def calculate_pad(self, *args, **kwargs):
        ...
    def calculate_rdf(self, *args, **kwargs):
        ...
    def calculate_vacf(self, traj: Trajectory, max_correlation_frames: typing.SupportsInt | typing.SupportsIndex = -1, start_frame: typing.SupportsInt | typing.SupportsIndex = 0, end_frame: typing.SupportsInt | typing.SupportsIndex = 18446744073709551615) -> None:
        """
        Calculate the Velocity Autocorrelation Function (VACF).
        
        Parameters
        ----------
        traj : Trajectory
            Trajectory with pre-calculated velocities.
        max_correlation_frames : int
            Max lag frames (-1 = half trajectory).
        start_frame : int
            First frame to include. Default 0.
        end_frame : int
            Last frame (exclusive). Default all.
        """
    def calculate_vdos(self) -> None:
        """
        Calculate the Vibrational Density of States (VDOS) from the VACF.
        Requires calculate_vacf() to have been called first.
        """
    def calculate_xrd(self, wavelength: typing.SupportsFloat | typing.SupportsIndex = 1.5406, theta_min: typing.SupportsFloat | typing.SupportsIndex = 5.0, theta_max: typing.SupportsFloat | typing.SupportsIndex = 90.0, bin_width: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
        Calculate the X-Ray Diffraction (XRD) pattern.
        
        Parameters
        ----------
        lambda : float
            X-ray wavelength (Å). Default 1.5406 (Cu Kα).
        theta_min : float
            Minimum 2θ angle (°). Default 5.0.
        theta_max : float
            Maximum 2θ angle (°). Default 90.0.
        bin_width : float
            2θ bin width (°). Default 1.0.
        """
    def get_all_histograms(self) -> dict[str, Histogram]:
        """
        Return a dict of all calculated histograms.
        """
    def get_ashcroft_weights(self) -> dict[str, float]:
        """
        Return the Ashcroft-Langreth weights used for S(Q) partials.
        """
    def get_available_histograms(self) -> list[str]:
        """
        Return a list of names for all currently available histograms.
        """
    def get_deborah_number(self) -> float:
        """
        Get the Deborah number.
        """
    def get_diffusion_coefficient_msd(self) -> float:
        """
        Get the self-diffusion coefficient computed from MSD (Å²/fs).
        """
    def get_diffusion_coefficient_vacf(self) -> float:
        """
        Get the self-diffusion coefficient computed from VACF (Å²/fs).
        """
    def get_histogram(self, name: str) -> Histogram:
        """
        Return the Histogram for the given key (e.g. 'g(r)', 'S(Q)').
        Raises KeyError if the histogram has not been calculated.
        """
    def get_relaxation_time(self) -> float:
        """
        Get the relaxation time computed from normalized VACF (fs).
        """
    def scale(self, factor: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Scale all histogram values by *factor* (used for normalization).
        """
    def smooth(self, name: str, sigma: typing.SupportsFloat | typing.SupportsIndex, kernel: KernelType = ...) -> None:
        """
        Smooth a specific histogram.
        
        Parameters
        ----------
        name : str
            Histogram name (e.g. 'g(r)').
        sigma : float
            Kernel bandwidth.
        kernel : KernelType
            Smoothing kernel. Default Gaussian.
        """
    def smooth_all(self, sigma: typing.SupportsFloat | typing.SupportsIndex, kernel: KernelType = ...) -> None:
        """
        Smooth all available histograms.
        
        sigma : float
            Kernel bandwidth.
        kernel : KernelType
            Smoothing kernel. Default Gaussian.
        """
class Element:
    """
    Chemical element specification.
    """
    def __init__(self) -> None:
        """
        Construct an empty Element.
        """
    @property
    def id(self) -> int:
        """
        Unique integer identifier for the element type.
        """
    @id.setter
    def id(self, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def symbol(self) -> str:
        """
        Chemical symbol (e.g., 'C', 'O', 'Fe').
        """
    @symbol.setter
    def symbol(self, arg0: str) -> None:
        ...
class GaussianRBFConfig:
    """
    Configuration parameters for Gaussian radial basis function expansion.
    """
    def __init__(self, start: typing.SupportsFloat | typing.SupportsIndex = 0.0, stop: typing.SupportsFloat | typing.SupportsIndex = 5.0, num_basis: typing.SupportsInt | typing.SupportsIndex = 8) -> None:
        ...
    @property
    def num_basis(self) -> int:
        """
        Number of Gaussian basis centers.
        """
    @num_basis.setter
    def num_basis(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def start(self) -> float:
        """
        Start center distance in Angstroms.
        """
    @start.setter
    def start(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def stop(self) -> float:
        """
        Stop center distance in Angstroms.
        """
    @stop.setter
    def stop(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class GraphDescriptors:
    """
    Extracts topological, structural, and spectral descriptors from PeriodicGraphData.
    """
    @staticmethod
    def compute_cna_descriptor(graph: PeriodicGraphData) -> list[int]:
        """
        Compute per-atom Common Neighbor Analysis (CNA) classification labels.
        """
    @staticmethod
    def compute_coordination_embedding(graph: PeriodicGraphData) -> list[float]:
        """
        Compute per-atom coordination number embedding.
        """
    @staticmethod
    def compute_graph_spectrum(graph: PeriodicGraphData, k: typing.SupportsInt | typing.SupportsIndex) -> list[float]:
        """
        Compute top-k eigenvalues of the graph adjacency matrix.
        """
    @staticmethod
    def compute_ring_statistics_descriptor(graph: PeriodicGraphData, max_size: typing.SupportsInt | typing.SupportsIndex = 6) -> list[float]:
        """
        Compute per-atom ring statistics embedding using cycle basis detection.
        """
    @staticmethod
    def populate_descriptors(graph: PeriodicGraphData, max_ring_size: typing.SupportsInt | typing.SupportsIndex = 6) -> None:
        """
        Populate all descriptor fields in PeriodicGraphData in-place.
        """
class Histogram:
    """
    Container for a single calculated distribution function.
    """
    def get_bins_numpy(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Return x-axis bins as a NumPy array (copy).
        """
    def get_partial_numpy(self, key: str) -> numpy.typing.NDArray[numpy.float32]:
        """
        Return a specific partial distribution as a NumPy array.
        """
    def get_smoothed_partial_numpy(self, key: str) -> numpy.typing.NDArray[numpy.float32]:
        """
        Return a specific smoothed partial distribution as a NumPy array.
        """
    @property
    def bins(self) -> list[float]:
        """
        List of bin center coordinates.
        """
    @bins.setter
    def bins(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def partials(self) -> dict[str, list[float]]:
        """
        Dict mapping partial key (e.g. 'Si-O') to y-values.
        """
    @partials.setter
    def partials(self, arg0: collections.abc.Mapping[str, collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]]) -> None:
        ...
    @property
    def smoothed_partials(self) -> dict[str, list[float]]:
        """
        Dict mapping partial key to smoothed y-values.
        """
    @smoothed_partials.setter
    def smoothed_partials(self, arg0: collections.abc.Mapping[str, collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]]) -> None:
        ...
    @property
    def title(self) -> str:
        """
        Descriptive title.
        """
    @property
    def x_label(self) -> str:
        """
        X-axis label.
        """
    @property
    def x_unit(self) -> str:
        """
        X-axis physical unit (e.g. 'Å').
        """
    @property
    def y_label(self) -> str:
        """
        Y-axis label.
        """
class KernelType:
    """
    Kernel type used for post-processing smoothing of histograms.
    
    Members:
    
      Gaussian : Gaussian (normal) kernel — smooth, infinite support.
    
      Bump : Infinitely-smooth bump function with compact support.
    
      Triweight : Triweight polynomial kernel with compact support.
    
      Epanechnikov : Optimal MISE kernel — parabolic, compact support.
    
      Cosine : Cosine kernel with compact support.
    
      Biweight : Quartic (biweight) kernel with compact support.
    """
    Biweight: typing.ClassVar[KernelType]  # value = <KernelType.Biweight: 5>
    Bump: typing.ClassVar[KernelType]  # value = <KernelType.Bump: 1>
    Cosine: typing.ClassVar[KernelType]  # value = <KernelType.Cosine: 4>
    Epanechnikov: typing.ClassVar[KernelType]  # value = <KernelType.Epanechnikov: 3>
    Gaussian: typing.ClassVar[KernelType]  # value = <KernelType.Gaussian: 0>
    Triweight: typing.ClassVar[KernelType]  # value = <KernelType.Triweight: 2>
    __members__: typing.ClassVar[dict[str, KernelType]]  # value = {'Gaussian': <KernelType.Gaussian: 0>, 'Bump': <KernelType.Bump: 1>, 'Triweight': <KernelType.Triweight: 2>, 'Epanechnikov': <KernelType.Epanechnikov: 3>, 'Cosine': <KernelType.Cosine: 4>, 'Biweight': <KernelType.Biweight: 5>}
    @typing.overload
    def __eq__(self, other: KernelType) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: KernelType) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class MLIPCalculator(BaseCalculator):
    """
    Calculator interface for machine learning interatomic potentials.
    """
    @staticmethod
    def calculate(cell: Cell, model: MLIPInterface | None = None) -> MLIPOutput:
        """
        Evaluate MLIP on a given cell.
        """
    def __init__(self) -> None:
        ...
class MLIPInterface:
    """
    Abstract interface for MLIP engines.
    """
    def __init__(self) -> None:
        ...
    def evaluate(self, cell: Cell) -> MLIPOutput:
        """
        Evaluate model on an atomic cell.
        """
    def get_model_name(self) -> str:
        """
        Return model descriptor name.
        """
class MLIPOutput:
    """
    Container for machine learning interatomic potential outputs.
    """
    def __init__(self) -> None:
        ...
    @property
    def forces(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Predicted atomic forces as a (N, 3) NumPy array.
        """
    @property
    def ldos(self) -> list[list[float]]:
        """
        Local Density of States matrix [N_atoms x N_bins].
        """
    @ldos.setter
    def ldos(self, arg0: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]]) -> None:
        ...
    @property
    def ldos_bins(self) -> int:
        """
        Number of LDOS energy bins.
        """
    @ldos_bins.setter
    def ldos_bins(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def per_atom_energy(self) -> list[float]:
        """
        Site-resolved per-atom energy.
        """
    @per_atom_energy.setter
    def per_atom_energy(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def total_energy(self) -> float:
        """
        Total predicted potential energy.
        """
    @total_energy.setter
    def total_energy(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class MotifProjectedTDOS:
    """
    Container for motif-partitioned Total Density of States spectra.
    """
    def __init__(self) -> None:
        ...
    def to_histogram(self, title: str = 'Motif-Projected TDOS') -> Histogram:
        """
        Convert motif-projected TDOS into a standard Histogram.
        """
    @property
    def energies(self) -> list[float]:
        """
        Energy grid values [eV].
        """
    @energies.setter
    def energies(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def frame_count(self) -> int:
        """
        Evaluated trajectory frame count.
        """
    @frame_count.setter
    def frame_count(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def motif_tdos(self) -> dict[str, list[float]]:
        """
        Partial TDOS spectra mapped by motif name.
        """
    @motif_tdos.setter
    def motif_tdos(self, arg0: collections.abc.Mapping[str, collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]]) -> None:
        ...
    @property
    def total_tdos(self) -> list[float]:
        """
        Total aggregated TDOS spectrum.
        """
    @total_tdos.setter
    def total_tdos(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
class OrbDescriptorConfig:
    """
    Configuration parameters for ORB-v3 graph and descriptor extraction.
    """
    def __init__(self, r_max: typing.SupportsFloat | typing.SupportsIndex = 6.0, num_rbf: typing.SupportsInt | typing.SupportsIndex = 8, l_max: typing.SupportsInt | typing.SupportsIndex = 3, include_self_loops: bool = False, compute_orb_features: bool = True) -> None:
        ...
    @property
    def compute_orb_features(self) -> bool:
        """
        Whether to compute fused outer-product [E, 128] ORB edge features.
        """
    @compute_orb_features.setter
    def compute_orb_features(self, arg0: bool) -> None:
        ...
    @property
    def include_self_loops(self) -> bool:
        """
        Whether to include self loops.
        """
    @include_self_loops.setter
    def include_self_loops(self, arg0: bool) -> None:
        ...
    @property
    def l_max(self) -> int:
        """
        Maximum spherical harmonics degree l.
        """
    @l_max.setter
    def l_max(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_rbf(self) -> int:
        """
        Number of Bessel radial basis functions.
        """
    @num_rbf.setter
    def num_rbf(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def r_max(self) -> float:
        """
        Radial cutoff distance in Angstroms.
        """
    @r_max.setter
    def r_max(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class PeriodicGraphBuilder:
    """
    Constructs periodic neighbor graphs for atomic GNN evaluation.
    """
    @staticmethod
    def build_graph(cell: Cell, cutoff_radius: typing.SupportsFloat | typing.SupportsIndex = 5.0, include_self_loops: bool = False, l_max: typing.SupportsInt | typing.SupportsIndex = 0) -> PeriodicGraphData:
        """
        Build periodic neighbor graph data for a unit cell.
        """
    @staticmethod
    def build_orb_graph(cell: Cell, config: OrbDescriptorConfig = ...) -> PeriodicGraphData:
        """
        Construct periodic neighbor graph and compute native ORB-v3 descriptors.
        """
    @staticmethod
    def compute_bessel_basis(distance: typing.SupportsFloat | typing.SupportsIndex, cutoff_radius: typing.SupportsFloat | typing.SupportsIndex, num_basis: typing.SupportsInt | typing.SupportsIndex) -> list[float]:
        """
        Compute spherical Bessel radial basis.
        """
    @staticmethod
    def compute_cutoff_envelope(distance: typing.SupportsFloat | typing.SupportsIndex, cutoff_radius: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
        Compute smooth polynomial cutoff envelope.
        """
    @staticmethod
    @typing.overload
    def compute_gaussian_rbf(distance: typing.SupportsFloat | typing.SupportsIndex, start: typing.SupportsFloat | typing.SupportsIndex, stop: typing.SupportsFloat | typing.SupportsIndex, num_basis: typing.SupportsInt | typing.SupportsIndex) -> list[float]:
        """
        Compute Gaussian radial basis functions.
        """
    @staticmethod
    @typing.overload
    def compute_gaussian_rbf(distance: typing.SupportsFloat | typing.SupportsIndex, config: GaussianRBFConfig) -> list[float]:
        """
        Compute Gaussian radial basis functions using configuration struct.
        """
    @staticmethod
    def compute_orb_bessel_basis(distance: typing.SupportsFloat | typing.SupportsIndex, config: OrbDescriptorConfig = ...) -> list[float]:
        """
        Compute ORB-v3 Bessel radial basis functions.
        """
    @staticmethod
    def compute_orb_cutoff_envelope(distance: typing.SupportsFloat | typing.SupportsIndex, cutoff_radius: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
        Compute ORB-v3 polynomial cutoff envelope of order p=4.
        """
    @staticmethod
    def compute_orb_spherical_harmonics(unit_vec: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], l_max: typing.SupportsInt | typing.SupportsIndex = 3) -> list[float]:
        """
        Compute ORB-v3 e3nn component-normalized spherical harmonics.
        """
    @staticmethod
    def get_atomic_number(symbol: str) -> int:
        """
        Return atomic number (Z) for element symbol.
        """
class PeriodicGraphData:
    """
    Container for periodic neighbor graph tensors ready for GNN model inference.
    """
    def __init__(self) -> None:
        ...
    @property
    def atom_count(self) -> int:
        """
        Total number of atoms N.
        """
    @property
    def atomic_numbers(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        Zero-copy access to atomic numbers (Z) as a (N,) NumPy array.
        """
    @property
    def cell(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to lattice matrix as a (3, 3) NumPy array.
        """
    @property
    def cna_labels(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        Zero-copy access to CNA classification labels as a (N,) NumPy array.
        """
    @property
    def coordination_desc(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to coordination numbers as a (N,) NumPy array.
        """
    @property
    def edge_count(self) -> int:
        """
        Total number of directed edges E.
        """
    @property
    def edge_cutoff_envelope(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to polynomial cutoff envelope values f_c(d) as a (E,) NumPy array.
        """
    @property
    def edge_distances(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to Euclidean edge distances ||r_ij|| as a (E,) NumPy array.
        """
    @property
    def edge_index(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        Zero-copy access to directed edge indices (COO format) as a (2, E) NumPy array.
        """
    @property
    def edge_orb_features(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to fused ORB-v3 edge features f_cut * (RBF (x) Y_lm) as a (E, num_rbf * num_sh) NumPy array.
        """
    @property
    def edge_radial_basis(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to edge radial basis features as a (E, num_rbf) NumPy array.
        """
    @property
    def edge_shifts(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to periodic cell integer shift vectors as a (E, 3) NumPy array.
        """
    @property
    def edge_spherical_harmonics(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to spherical harmonics features as a (E, (l_max+1)^2) NumPy array.
        """
    @property
    def edge_unit_vectors(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to normalized unit displacement vectors r_hat_ij as a (E, 3) NumPy array.
        """
    @property
    def edge_vectors(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to Cartesian displacement vectors r_ij as a (E, 3) NumPy array.
        """
    @property
    def positions(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to atomic Cartesian coordinates as a (N, 3) NumPy array.
        """
    @property
    def ring_desc(self) -> numpy.typing.NDArray[numpy.float32]:
        """
        Zero-copy access to per-atom ring descriptors as a (N, max_ring_size) NumPy array.
        """
class StructuralElectronicCorrelation:
    """
    Trajectory-level structural-electronic correlation pipeline.
    """
    @staticmethod
    def correlate_cna(dists: DistributionFunctions, traj: Trajectory, params: TDOSParams) -> MotifProjectedTDOS:
        """
        Correlate trajectory LDoS with Common Neighbor Analysis (CNA) classifications.
        """
    @staticmethod
    def correlate_steinhardt(dists: DistributionFunctions, traj: Trajectory, params: TDOSParams) -> MotifProjectedTDOS:
        """
        Correlate trajectory LDoS with Steinhardt bond-order parameters.
        """
class StructureAnalyzer:
    """
    Computes pairwise distances, bond angles, and dihedral angles
    for a single simulation cell (frame).
    
    The tensors are indexed by element type: distances[e1][e2][pair_idx],
    angles[center][e1][e2][angle_idx], dihedrals[e1][e2][e3][e4][idx].
    """
    def __init__(self, cell: Cell, cutoff: typing.SupportsFloat | typing.SupportsIndex, bond_cutoffs_sq: typing.Any, ignore_periodic_self_interactions: bool = True) -> None:
        """
        Construct and immediately compute all pair data.
        
        Parameters
        ----------
        cell : Cell
            The periodic simulation cell.
        cutoff : float
            Neighbor search cutoff radius (Å).
        bond_cutoffs_sq : list[list[BondCutoffRange]] | list[list[tuple]] | list[list[float]]
            Per-element-pair bond cutoffs.
        ignore_periodic_self_interactions : bool
            If True, atoms do not interact with their own periodic images.
        """
    def angles(self) -> list[list[list[list[float]]]]:
        """
        4D tensor of bond angles [center_e][e1][e2][angle_idx].
        """
    def dihedrals(self) -> list[list[list[list[list[float]]]]]:
        """
        5D tensor of dihedral angles [e1][e2][e3][e4][dihedral_idx].
        """
    def raw_histograms(self) -> list[list[list[float]]]:
        """
        3D tensor of raw distance histogram bins [e1][e2][bin_idx].
        """
class TDOSCalculator(BaseCalculator):
    """
    Calculator for MLIP Total Density of States (TDOS).
    """
    @staticmethod
    def calculate(cell: Cell, params: TDOSParams = ...) -> Histogram:
        """
        Calculate Total Density of States for a single cell.
        """
    @staticmethod
    def calculate_trajectory(traj: Trajectory, params: TDOSParams = ...) -> Histogram:
        """
        Calculate frame-averaged Total Density of States for a trajectory.
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, model: MLIPInterface | None = None) -> None:
        """
        Construct TDOSCalculator with an optional MLIP model engine.
        """
    def get_model(self) -> MLIPInterface:
        """
        Get pointer to attached MLIP model engine.
        """
    def set_model(self, model: MLIPInterface) -> None:
        """
        Set or update the MLIP model engine.
        """
class TDOSParams:
    """
    Configuration parameters for Total Density of States (TDOS) calculation.
    """
    def __init__(self, e_min: typing.SupportsFloat | typing.SupportsIndex = -15.0, e_max: typing.SupportsFloat | typing.SupportsIndex = 5.0, model: MLIPInterface | None = None) -> None:
        ...
    @property
    def e_max(self) -> float:
        """
        Upper energy bound in eV (relative to E_F).
        """
    @e_max.setter
    def e_max(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def e_min(self) -> float:
        """
        Lower energy bound in eV (relative to E_F).
        """
    @e_min.setter
    def e_min(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def model(self) -> MLIPInterface:
        """
        MLIPInterface model pointer.
        """
    @model.setter
    def model(self, arg0: MLIPInterface) -> None:
        ...
class Trajectory:
    """
    Time-series collection of Cell simulation snapshots.
    """
    @staticmethod
    def from_ase(images: typing.Iterable[typing.Any]) -> typing.Any:
        """
        
        Convert an iterable or sequence of ASE Atoms frames into a Correlation Trajectory.
        
        Parameters
        ----------
        images : Iterable[ase.Atoms]
            Sequence of atomic snapshots (e.g., from ase.io.Trajectory).
        
        Returns
        -------
        correlation.Trajectory
            Constructed Trajectory containing all frames.
        """
    def to_ase(self) -> list[typing.Any]:
        """
        
        Convert a Correlation Trajectory into a list of ASE Atoms objects.
        
        Parameters
        ----------
        traj : correlation.Trajectory
            Source multi-frame trajectory.
        
        Returns
        -------
        list[ase.Atoms]
            List of converted ASE Atoms snapshots.
        """
    def __getitem__(self, index: typing.SupportsInt | typing.SupportsIndex) -> Cell:
        """
        Access a frame snapshot by integer index.
        """
    def __init__(self) -> None:
        """
        Construct an empty trajectory.
        """
    def __iter__(self) -> collections.abc.Iterator[Cell]:
        """
        Iterate over frames in the trajectory.
        """
    def __len__(self) -> int:
        """
        Number of frames in the trajectory.
        """
    def add_frame(self, frame: Cell) -> None:
        """
        Append a Cell frame to the trajectory.
        """
    def append(self, frame: Cell) -> None:
        """
        Append a Cell frame to the trajectory (alias for add_frame).
        """
    def num_frames(self) -> int:
        """
        Number of frames in the trajectory.
        """
    @property
    def frames(self) -> list[Cell]:
        """
        Direct reference to the sequence of Cell frames.
        """
    @property
    def time_step(self) -> float:
        """
        Time interval between consecutive frames.
        """
    @time_step.setter
    def time_step(self, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class TrajectoryAnalyzer:
    """
    Orchestrates structural analysis across multiple frames of a trajectory.
    
    Provides per-frame StructureAnalyzer factories and trajectory metadata.
    """
    def __init__(self, trajectory: Trajectory, neighbor_cutoff: typing.SupportsFloat | typing.SupportsIndex, bond_cutoffs: typing.Any, start_frame: typing.SupportsInt | typing.SupportsIndex = 0, end_frame: typing.SupportsInt | typing.SupportsIndex = -1, ignore_periodic_self_interactions: bool = True, progress_callback: collections.abc.Callable[[float, str], None] | None = None) -> None:
        """
        Construct a TrajectoryAnalyzer.
        
        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to analyze.
        neighbor_cutoff : float
            Global neighbor search cutoff radius (Å).
        bond_cutoffs : list[list[BondCutoffRange]] | list[list[tuple]] | list[list[float]]
            Per-element-pair bond cutoffs (Å).
        start_frame : int, optional
            Index of the first frame to analyze. Default 0.
        end_frame : int, optional
            Index of the last frame (-1 for all). Default -1.
        ignore_periodic_self_interactions : bool, optional
            Skip atom–own-image interactions. Default True.
        progress_callback : callable, optional
            Called with (fraction: float, message: str) during computation.
        """
    def create_analyzer(self, frame_idx: typing.SupportsInt | typing.SupportsIndex) -> StructureAnalyzer:
        """
        Create a StructureAnalyzer for the given frame index.
        """
    def get_neighbor_cutoff(self) -> float:
        """
        Global neighbor search cutoff radius (Å).
        """
    def get_num_frames(self) -> int:
        """
        Total number of frames in the analysis window.
        """
    def get_start_frame(self) -> int:
        """
        Index of the first frame being analyzed.
        """
    def get_time_step(self) -> float:
        """
        Time step between frames (from trajectory metadata).
        """
def build_orb_graph(cell: Cell, config: OrbDescriptorConfig = ...) -> PeriodicGraphData:
    """
    Convenience helper to construct periodic neighbor graph and extract ORB-v3 descriptors.
    """
def build_periodic_graph(cell: Cell, cutoff_radius: typing.SupportsFloat | typing.SupportsIndex = 5.0, include_self_loops: bool = False, l_max: typing.SupportsInt | typing.SupportsIndex = 0) -> PeriodicGraphData:
    """
    Convenience helper to construct periodic neighbor graph for GNN evaluation.
    """
def compute_cna_descriptor(graph: PeriodicGraphData) -> list[int]:
    """
    Compute per-atom Common Neighbor Analysis (CNA) classification labels.
    """
def compute_coordination_embedding(graph: PeriodicGraphData) -> list[float]:
    """
    Compute per-atom coordination number embedding.
    """
def compute_graph_spectrum(graph: PeriodicGraphData, k: typing.SupportsInt | typing.SupportsIndex) -> list[float]:
    """
    Compute top-k eigenvalues of the graph adjacency matrix.
    """
def compute_ring_statistics_descriptor(graph: PeriodicGraphData, max_size: typing.SupportsInt | typing.SupportsIndex = 6) -> list[float]:
    """
    Compute per-atom ring statistics embedding using cycle basis detection.
    """
def correlate_cna(dists: DistributionFunctions, traj: Trajectory, params: TDOSParams) -> MotifProjectedTDOS:
    """
    Correlate trajectory LDoS with Common Neighbor Analysis (CNA) classifications.
    """
def correlate_steinhardt(dists: DistributionFunctions, traj: Trajectory, params: TDOSParams) -> MotifProjectedTDOS:
    """
    Correlate trajectory LDoS with Steinhardt bond-order parameters.
    """
def get_all_calculators() -> list[BaseCalculator]:
    """
    Return a list of all registered calculator objects.
    """
def get_calculator(name: str) -> BaseCalculator:
    """
    Look up a registered calculator by its short name (e.g. 'RDF').
    
    Returns a reference into the factory's singleton registry.
    
    Parameters
    ----------
    name : str
        Short name of the desired calculator.
    
    Raises
    ------
    RuntimeError
        If no calculator with that name is registered.
    """
def get_writer(name: str) -> BaseWriter:
    """
    Look up a file writer by name (e.g. 'CSV', 'HDF5').
    
    Returns a BaseWriter reference into the factory's singleton registry.
    
    Parameters
    ----------
    name : str
        Format name of the desired writer.
    
    Raises
    ------
    RuntimeError
        If no writer with that name is registered.
    """
def get_writer_for_extension(extension: str) -> BaseWriter:
    """
    Look up a writer by file extension (e.g. '.csv').
    """
def list_calculators() -> list[str]:
    """
    Return a list of short names for all registered calculators
    (e.g. ['RDF', 'SQ', 'PAD', ...]).
    """
def list_writers() -> list[str]:
    """
    Return the names of all registered file writers.
    """
def populate_descriptors(graph: PeriodicGraphData, max_ring_size: typing.SupportsInt | typing.SupportsIndex = 6) -> None:
    """
    Populate all descriptor fields in PeriodicGraphData in-place.
    """
def read(filepath: str) -> Trajectory:
    """
    Read an atomic trajectory or structure file into a Trajectory object.
    
    Automatically detects format from file extension (e.g. .car, .arc, .lammps, .xyz).
    
    Parameters
    ----------
    filepath : str
        Path to the input structure or trajectory file.
    
    Returns
    -------
    Trajectory
        Parsed trajectory containing one or more simulation Cell frames.
    
    Raises
    ------
    RuntimeError
        If no reader is available for the given file extension or parsing fails.
    """
def write_csv(base_path: str, dists: DistributionFunctions, write_smoothed: bool = False) -> None:
    """
    Convenience function: write all histograms in *dists* to CSV files.
    
    Parameters
    ----------
    base_path : str
        Base name for output files (e.g. 'output/sample').
    dists : DistributionFunctions
        The analysis results to export.
    write_smoothed : bool, optional
        If True, also write smoothed data files. Default is False.
    """
BCC: CNALabel  # value = <CNALabel.BCC: 3>
Biweight: KernelType  # value = <KernelType.Biweight: 5>
Bump: KernelType  # value = <KernelType.Bump: 1>
Cosine: KernelType  # value = <KernelType.Cosine: 4>
Epanechnikov: KernelType  # value = <KernelType.Epanechnikov: 3>
FCC: CNALabel  # value = <CNALabel.FCC: 1>
Gaussian: KernelType  # value = <KernelType.Gaussian: 0>
HCP: CNALabel  # value = <CNALabel.HCP: 2>
ICO: CNALabel  # value = <CNALabel.ICO: 4>
OTHER: CNALabel  # value = <CNALabel.OTHER: 0>
Triweight: KernelType  # value = <KernelType.Triweight: 2>
