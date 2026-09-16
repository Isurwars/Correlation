"""

Correlation — High-performance structural analysis for atomistic simulations.

This package provides Python bindings for the Correlation C++ analysis engine,
enabling computation of radial distribution functions (RDF), structure factors
S(Q), bond angle distributions, and other structural properties of liquids,
amorphous solids, and crystalline materials.

Example usage::

    import correlation

    # Read a structure file
    cell = correlation.Cell.from_file("structure.car")

    # Compute distribution functions
    df = correlation.DistributionFunctions(cell, cutoff=10.0)
    df.calculate_rdf(r_max=20.0, bin_width=0.05)

    # Access results
    hist = df.get_histogram("g_r")
    print(hist.bins, hist.partials)
"""
from __future__ import annotations
from correlation._correlation import AnalysisSettings
from correlation._correlation import Atom as _orig_atom
from correlation._correlation import BaseCalculator
from correlation._correlation import BaseWriter
from correlation._correlation import BondCutoffRange
from correlation._correlation import CNALabel
from correlation._correlation import CSVWriter
from correlation._correlation import Cell
from correlation._correlation import DistributionFunctions
from correlation._correlation import Element
from correlation._correlation import GaussianRBFConfig
from correlation._correlation import GraphDescriptors
from correlation._correlation import Histogram
from correlation._correlation import KernelType
from correlation._correlation import MLIPCalculator
from correlation._correlation import MLIPInterface
from correlation._correlation import MLIPOutput
from correlation._correlation import MotifProjectedTDOS
from correlation._correlation import OrbDescriptorConfig
from correlation._correlation import PeriodicGraphBuilder
from correlation._correlation import PeriodicGraphData
from correlation._correlation import StructuralElectronicCorrelation
from correlation._correlation import StructureAnalyzer
from correlation._correlation import TDOSCalculator
from correlation._correlation import TDOSParams
from correlation._correlation import Trajectory
from correlation._correlation import TrajectoryAnalyzer
from correlation._correlation import build_orb_graph
from correlation._correlation import build_periodic_graph
from correlation._correlation import compute_cna_descriptor
from correlation._correlation import compute_coordination_embedding
from correlation._correlation import compute_graph_spectrum
from correlation._correlation import compute_ring_statistics_descriptor
from correlation._correlation import correlate_cna
from correlation._correlation import correlate_steinhardt
from correlation._correlation import get_all_calculators
from correlation._correlation import get_calculator
from correlation._correlation import get_writer
from correlation._correlation import get_writer_for_extension
from correlation._correlation import list_calculators as get_registered_calculators
from correlation._correlation import list_calculators
from correlation._correlation import list_writers
from correlation._correlation import populate_descriptors
from correlation._correlation import read
from correlation._correlation import write_csv
from correlation.adapters import from_ase
from correlation.adapters import from_ase_trajectory
from correlation.adapters import from_pymatgen
from correlation.adapters import to_ase
from correlation.adapters import to_ase_trajectory
from correlation.adapters import to_pymatgen
from . import _correlation
from . import adapters
__all__: list[str] = ['AnalysisSettings', 'Atom', 'AtomShim', 'BCC', 'BaseCalculator', 'BaseWriter', 'Biweight', 'BondCutoffRange', 'Bump', 'CNALabel', 'CSVWriter', 'Cell', 'Cosine', 'DistributionFunctions', 'Element', 'Epanechnikov', 'FCC', 'Gaussian', 'GaussianRBFConfig', 'GraphDescriptors', 'HCP', 'Histogram', 'ICO', 'KernelType', 'MLIPCalculator', 'MLIPInterface', 'MLIPOutput', 'MotifProjectedTDOS', 'OTHER', 'OrbDescriptorConfig', 'PeriodicGraphBuilder', 'PeriodicGraphData', 'RDFParams', 'StructuralElectronicCorrelation', 'StructureAnalyzer', 'TDOSCalculator', 'TDOSParams', 'Trajectory', 'TrajectoryAnalyzer', 'Triweight', 'adapters', 'build_orb_graph', 'build_periodic_graph', 'compute_cna_descriptor', 'compute_coordination_embedding', 'compute_graph_spectrum', 'compute_ring_statistics_descriptor', 'correlate_cna', 'correlate_steinhardt', 'from_ase', 'from_ase_trajectory', 'from_pymatgen', 'get_all_calculators', 'get_calculator', 'get_registered_calculators', 'get_writer', 'get_writer_for_extension', 'list_calculators', 'list_writers', 'populate_descriptors', 'read', 'to_ase', 'to_ase_trajectory', 'to_pyg', 'to_pymatgen', 'to_torch_geometric', 'write_csv']
class AtomShim:
    """
    Ergonomic container for atom specifications.
    """
    def __init__(self, symbol: str, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        ...
    def __repr__(self) -> str:
        ...
    @property
    def position(self) -> list[float]:
        ...
    @property
    def symbol(self) -> str:
        ...
class RDFParams:
    """
    Parameters container for radial distribution function analysis.
    """
    def __init__(self, r_max: float = 20.0, r_bin_width: float = 0.05, **kwargs):
        ...
    def __repr__(self) -> str:
        ...
class _CallableInt(int):
    """
    Integer that can also be invoked as a zero-argument callable.
    """
    def __call__(self) -> int:
        ...
def _atom_constructor(*args, **kwargs):
    ...
def _cell_add_atom(self, *args, **kwargs):
    ...
def _cell_init(self, *args, **kwargs):
    ...
def _df_calc_pad(self, *args, **kwargs):
    ...
def _df_calc_rdf(self, *args, **kwargs):
    ...
def _df_init(self, cell, cutoff: float = 0.0, bond_cutoffs = None, **kwargs):
    ...
def to_torch_geometric(graph_data):
    """
    
    Converts a PeriodicGraphData instance into a PyTorch Geometric Data object.
    
    Parameters
    ----------
    graph_data : PeriodicGraphData
        Periodic neighbor graph constructed by build_periodic_graph() or build_orb_graph().
    
    Returns
    -------
    torch_geometric.data.Data
        PyG graph data containing tensors:
        - z: Atomic numbers (N,)
        - pos: Atomic Cartesian coordinates (N, 3)
        - edge_index: Graph directed edges (2, E)
        - edge_shift: Periodic shift vectors (E, 3)
        - edge_vec: Cartesian displacement vectors r_ij (E, 3)
        - edge_dist: Euclidean distances ||r_ij|| (E,)
        - cell: Lattice vector matrix (3, 3)
        - pbc: Periodic boundary flags [True, True, True]
        - edge_unit_vec: Normalized unit displacement vectors r_hat (E, 3) [if present]
        - edge_rbf: Bessel radial basis expansion (E, 8) [if present]
        - edge_sh: e3nn spherical harmonics (E, 16) [if present]
        - edge_cutoff: Order p=4 polynomial cutoff envelope (E,) [if present]
        - edge_orb_features: Fused ORB-v3 edge attributes (E, 128) [if present]
    """
BCC: _correlation.CNALabel  # value = <CNALabel.BCC: 3>
Biweight: _correlation.KernelType  # value = <KernelType.Biweight: 5>
Bump: _correlation.KernelType  # value = <KernelType.Bump: 1>
Cosine: _correlation.KernelType  # value = <KernelType.Cosine: 4>
Epanechnikov: _correlation.KernelType  # value = <KernelType.Epanechnikov: 3>
FCC: _correlation.CNALabel  # value = <CNALabel.FCC: 1>
Gaussian: _correlation.KernelType  # value = <KernelType.Gaussian: 0>
HCP: _correlation.CNALabel  # value = <CNALabel.HCP: 2>
ICO: _correlation.CNALabel  # value = <CNALabel.ICO: 4>
OTHER: _correlation.CNALabel  # value = <CNALabel.OTHER: 0>
Triweight: _correlation.KernelType  # value = <KernelType.Triweight: 2>
_DF_CELL_MAP: dict = {}
_orig_atom_count: property  # value = <property object>
Atom = _atom_constructor
to_pyg = to_torch_geometric
