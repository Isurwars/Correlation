"""
test_disorder_descriptors.py — Validation of Python bindings for disorder and order parameter descriptors.

Tests:
  - Ten Wolde-Frenkel (TWF) mutual bond coherence (N_conn)
  - Chau-Hardwick tetrahedral order parameter (q_tet)
  - Local Entropic Fingerprints (s2)
  - Radial bond correlation g6(r) & correlation length fitting
  - Integration with DistributionFunctions
"""

import math
import sys
import numpy as np
import pytest

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import correlation
else:
    try:
        import correlation
    except ImportError:
        try:
            import matcorr as correlation
        except ImportError:
            correlation = None

if correlation is None or not hasattr(correlation, "Cell"):
    pytest.skip(
        "correlation/matcorr C++ extension module not found (run 'pip install .' to test bindings)",
        allow_module_level=True,
    )


def create_fcc_cu_cell():
    """Create a 2x2x2 FCC supercell of Cu (32 atoms, a = 3.615 A)."""
    a = 3.615
    cell = correlation.Cell(
        [2.0 * a, 0.0, 0.0],
        [0.0, 2.0 * a, 0.0],
        [0.0, 0.0, 2.0 * a],
    )
    # Basis vectors of FCC
    basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ]
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                for bx, by, bz in basis:
                    pos = [
                        (ix + bx) * a,
                        (iy + by) * a,
                        (iz + bz) * a,
                    ]
                    cell.add_atom("Cu", pos)
    return cell


def create_tetrahedral_cell():
    """Create a non-periodic box containing a central atom surrounded tetrahedrally by 4 neighbors."""
    cell = correlation.Cell(
        [20.0, 0.0, 0.0],
        [0.0, 20.0, 0.0],
        [0.0, 0.0, 20.0],
    )
    center = [10.0, 10.0, 10.0]
    cell.add_atom("Si", center)
    # Standard tetrahedron vertices scaled to bond length ~1.6 A
    d = 1.6 / math.sqrt(3.0)
    offsets = [
        [ d,  d,  d],
        [ d, -d, -d],
        [-d,  d, -d],
        [-d, -d,  d],
    ]
    for dx, dy, dz in offsets:
        cell.add_atom("O", [center[0] + dx, center[1] + dy, center[2] + dz])
    return cell


def test_ten_wolde_frenkel_fcc():
    """Verify TWF N_conn detects high coherence in ideal crystalline FCC Cu."""
    cell = create_fcc_cu_cell()
    assert len(cell.atoms) == 32

    analyzer = correlation.StructureAnalyzer(cell, 3.0, [[9.0]], True)
    params = correlation.TenWoldeFrenkelParams(0.7, 7)
    assert params.chi_cut == pytest.approx(0.7)
    assert params.nconn_min == 7

    n_conn = correlation.compute_ten_wolde_frenkel(cell, analyzer, params)
    assert isinstance(n_conn, np.ndarray)
    assert n_conn.shape == (32,)
    # In an ideal FCC lattice, all 12 nearest neighbors share identical orientation (chi_ij ~ 1.0)
    for val in n_conn:
        assert val >= 7.0
        assert val <= 12.0


def test_tetrahedral_order_ideal():
    """Verify q_tet gives 1.0 for an ideal tetrahedral configuration."""
    cell = create_tetrahedral_cell()
    assert len(cell.atoms) == 5

    analyzer = correlation.StructureAnalyzer(cell, 2.0, [[4.0, 4.0], [4.0, 4.0]], True)
    params = correlation.TetrahedralParams()

    q_tet = correlation.compute_tetrahedral_order(cell, analyzer, params)
    assert isinstance(q_tet, np.ndarray)
    assert q_tet.shape == (5,)
    # The central Si atom (index 0) has 4 bonded O neighbors with 109.47 degree angles -> q_tet = 1.0
    assert q_tet[0] == pytest.approx(1.0, abs=1e-3)


def test_local_entropy_fcc():
    """Verify Local Entropy (s2) produces negative entropy values for dense crystal packing."""
    cell = create_fcc_cu_cell()
    analyzer = correlation.StructureAnalyzer(cell, 3.0, [[9.0]], True)

    params = correlation.LocalEntropyParams(cutoff=5.0, sigma=0.2)
    assert params.cutoff == pytest.approx(5.0)
    assert params.sigma == pytest.approx(0.2)

    s2 = correlation.compute_local_entropy(cell, analyzer, params)
    assert isinstance(s2, np.ndarray)
    assert s2.shape == (32,)
    assert np.all(np.isfinite(s2))
    # Two-body excess entropy is strictly negative for ordered crystals
    assert np.all(s2 < 0.0)


def test_g6_and_correlation_fit():
    """Verify g6(r) calculation and log-linear correlation length exponential fitting."""
    # 1. Exponential decay fit verification
    target_xi = 3.5
    r_bins = np.linspace(1.0, 10.0, 20)
    g6_vals = np.exp(-r_bins / target_xi)

    fit = correlation.fit_correlation_length(r_bins, g6_vals)
    assert not fit.is_plateau
    assert fit.xi_6 == pytest.approx(target_xi, abs=0.05)
    assert fit.r_squared > 0.99

    # 2. FCC supercell g6(r) evaluation
    cell = create_fcc_cu_cell()
    analyzer = correlation.StructureAnalyzer(cell, 3.0, [[9.0]], True)
    hist = correlation.compute_g6(cell, analyzer, r_max=5.0, r_bin_width=0.1)

    assert hist.title == "g₆(r) — Bond-Orientational Correlation Function"
    assert "Total" in hist.partials
    assert "Counts" in hist.partials
    assert len(hist.bins) == len(hist.partials["Total"])


def test_distribution_functions_integration():
    """Verify all 4 new calculators can be executed via DistributionFunctions."""
    cell = create_fcc_cu_cell()
    df = correlation.DistributionFunctions(cell, cutoff=3.0, bond_cutoffs=[[9.0]])

    df.calculate_twf(chi_cut=0.7, nconn_min=7)
    df.calculate_tetrahedral_order()
    df.calculate_local_entropy(cutoff=5.0, sigma=0.2)
    df.calculate_g6(r_max=5.0, r_bin_width=0.1)

    histograms = df.get_all_histograms()
    assert "TWF_Nconn" in histograms
    assert "Q_tet" in histograms
    assert "s2" in histograms
    assert "G6_r" in histograms
