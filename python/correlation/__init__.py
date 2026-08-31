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
        Periodic neighbor graph constructed by build_periodic_graph().

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
    """
    try:
        import torch
        from torch_geometric.data import Data
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

    return Data(
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


to_pyg = to_torch_geometric

