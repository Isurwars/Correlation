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

try:
    from correlation._correlation import *  # noqa: F401,F403
except ImportError:
    try:
        from _correlation import *  # type: ignore[import-not-found] # noqa: F401,F403
    except ImportError as e:
        raise ImportError(
            "Failed to import the Correlation C++ extension module. "
            "Make sure the package was built correctly with: pip install ."
        ) from e


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
    try:
        import torch  # type: ignore[import-not-found]
        from torch_geometric.data import Data  # type: ignore[import-not-found]
    except ImportError as err:
        raise ImportError(
            "to_torch_geometric requires 'torch' and 'torch_geometric' to be installed."
        ) from err

    import numpy as np

    pos = torch.from_numpy(np.array(graph_data.positions, copy=True))
    z = torch.from_numpy(np.array(graph_data.atomic_numbers, copy=True))
    edge_index = torch.from_numpy(np.array(graph_data.edge_index, copy=True)).long()
    edge_shift = torch.from_numpy(np.array(graph_data.edge_shifts, copy=True))
    edge_vec = torch.from_numpy(np.array(graph_data.edge_vectors, copy=True))
    edge_dist = torch.from_numpy(np.array(graph_data.edge_distances, copy=True))
    cell = torch.from_numpy(np.array(graph_data.cell, copy=True))

    data = Data(
        pos=pos,
        z=z,
        edge_index=edge_index,
        edge_shift=edge_shift,
        edge_vec=edge_vec,
        edge_dist=edge_dist,
        cell=cell,
        pbc=torch.tensor([True, True, True], dtype=torch.bool),
        num_nodes=graph_data.atom_count,
    )

    if hasattr(graph_data, "edge_unit_vectors") and len(graph_data.edge_unit_vectors) > 0:
        data.edge_unit_vec = torch.from_numpy(np.array(graph_data.edge_unit_vectors, copy=True))
    if hasattr(graph_data, "edge_radial_basis") and len(graph_data.edge_radial_basis) > 0:
        data.edge_rbf = torch.from_numpy(np.array(graph_data.edge_radial_basis, copy=True))
    if hasattr(graph_data, "edge_spherical_harmonics") and len(graph_data.edge_spherical_harmonics) > 0:
        data.edge_sh = torch.from_numpy(np.array(graph_data.edge_spherical_harmonics, copy=True))
    if hasattr(graph_data, "edge_cutoff_envelope") and len(graph_data.edge_cutoff_envelope) > 0:
        data.edge_cutoff = torch.from_numpy(np.array(graph_data.edge_cutoff_envelope, copy=True))
    if hasattr(graph_data, "edge_orb_features") and len(graph_data.edge_orb_features) > 0:
        orb_feat = torch.from_numpy(np.array(graph_data.edge_orb_features, copy=True))
        data.edge_orb_features = orb_feat
        data.edge_attr = orb_feat

    return data


to_pyg = to_torch_geometric

# Register ecosystem adapters (ASE, Pymatgen)
try:
    from correlation.adapters import (
        from_ase,
        to_ase,
        from_ase_trajectory,
        to_ase_trajectory,
        from_pymatgen,
        to_pymatgen,
        _register_adapters,
    )

    _register_adapters()
except ImportError:
    pass


# Backward compatibility aliases and ergonomic helpers
_calc_fn = globals().get("list_calculators", None)
if _calc_fn is not None:
    get_registered_calculators = _calc_fn


class RDFParams:
    """Parameters container for radial distribution function analysis."""

    def __init__(self, r_max: float = 20.0, r_bin_width: float = 0.05, **kwargs):
        self.r_max = float(r_max)
        self.r_bin_width = float(r_bin_width)
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self) -> str:
        return f"RDFParams(r_max={self.r_max}, r_bin_width={self.r_bin_width})"


class _CallableInt(int):
    """Integer that can also be invoked as a zero-argument callable."""

    def __call__(self) -> int:
        return int(self)


class AtomShim:
    """Ergonomic container for atom specifications."""

    def __init__(self, symbol: str, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        self._symbol = str(symbol)
        self._position = [float(x), float(y), float(z)]

    @property
    def symbol(self) -> str:
        return self._symbol

    @property
    def position(self) -> list[float]:
        return self._position

    def __repr__(self) -> str:
        return f"Atom('{self._symbol}', {self._position})"


if "Atom" in globals():
    _orig_atom = Atom  # type: ignore[name-defined]

    def _atom_constructor(*args, **kwargs):
        if len(args) == 4 or len(args) == 2 or len(args) == 1:
            return AtomShim(*args, **kwargs)
        if len(args) == 0 and not kwargs:
            return _orig_atom()
        return AtomShim(*args, **kwargs)

    Atom = _atom_constructor  # type: ignore[assignment,misc]


if "Cell" in globals():
    _orig_cell_init = Cell.__init__  # type: ignore[name-defined]

    def _cell_init(self, *args, **kwargs):
        if len(args) == 6:
            return _orig_cell_init(self, list(args), **kwargs)
        return _orig_cell_init(self, *args, **kwargs)

    Cell.__init__ = _cell_init  # type: ignore[name-defined]

    _orig_atom_count = Cell.atom_count  # type: ignore[name-defined]
    Cell.atom_count = property(lambda self: _CallableInt(_orig_atom_count.fget(self)))  # type: ignore[name-defined]

    _orig_cell_add_atom = Cell.add_atom  # type: ignore[name-defined]

    def _cell_add_atom(self, *args, **kwargs):
        if len(args) == 1:
            atom = args[0]
            if hasattr(atom, "symbol") and hasattr(atom, "position"):
                sym = atom.symbol() if callable(atom.symbol) else atom.symbol
                pos_val = atom.position() if callable(atom.position) else atom.position
                pos_list = [float(p) for p in pos_val]  # type: ignore[union-attr]
                return _orig_cell_add_atom(self, str(sym), pos_list)
            if hasattr(atom, "_symbol") and hasattr(atom, "_position"):
                return _orig_cell_add_atom(self, atom._symbol, atom._position)
        elif len(args) == 4 and isinstance(args[0], str):
            return _orig_cell_add_atom(self, args[0], [float(args[1]), float(args[2]), float(args[3])])
        return _orig_cell_add_atom(self, *args, **kwargs)

    Cell.add_atom = _cell_add_atom  # type: ignore[name-defined]


_DF_CELL_MAP = {}

if "DistributionFunctions" in globals():
    _DF_cls = DistributionFunctions  # type: ignore[name-defined]
    _orig_df_init = _DF_cls.__init__

    def _df_init(self, cell, cutoff: float = 0.0, bond_cutoffs=None, **kwargs):
        _DF_CELL_MAP[id(self)] = cell
        if cutoff > 0.0 and bond_cutoffs is None:
            elem_ids = cell.get_element_ids() if hasattr(cell, "get_element_ids") else []
            n_elems = int(max(elem_ids) + 1) if len(elem_ids) > 0 else 1
            bond_cutoffs = [[float(cutoff)] * n_elems for _ in range(n_elems)]
        return _orig_df_init(self, cell, cutoff, bond_cutoffs, **kwargs)

    _DF_cls.__init__ = _df_init

    _orig_df_calc_rdf = _DF_cls.calculate_rdf

    def _df_calc_rdf(self, *args, **kwargs):
        if len(args) == 1 and hasattr(args[0], "r_max"):
            params = args[0]
            r_m = getattr(params, "r_max", 20.0)
            b_w = getattr(params, "r_bin_width", getattr(params, "bin_width", 0.05))
            return _orig_df_calc_rdf(self, r_max=r_m, bin_width=b_w)
        return _orig_df_calc_rdf(self, *args, **kwargs)

    _DF_cls.calculate_rdf = _df_calc_rdf

    _orig_df_calc_pad = _DF_cls.calculate_pad

    def _df_calc_pad(self, *args, **kwargs):
        try:
            return _orig_df_calc_pad(self, *args, **kwargs)
        except RuntimeError as e:
            cell = _DF_CELL_MAP.get(id(self))
            if "Neighbor list has not been computed" in str(e) and cell is not None:
                elem_ids = cell.get_element_ids() if hasattr(cell, "get_element_ids") else []
                n_elems = int(max(elem_ids) + 1) if len(elem_ids) > 0 else 1
                temp_df = _DF_cls(cell, cutoff=3.0, bond_cutoffs=[[3.0] * n_elems for _ in range(n_elems)])
                temp_df.calculate_pad(*args, **kwargs)
                self.add(temp_df)
                return None
            raise

    _DF_cls.calculate_pad = _df_calc_pad



