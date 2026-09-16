"""
test.py — End-to-end smoke test for the correlation Python bindings.

Run from the build directory after building the correlation_py target:
    cmake --build . --target correlation_py -j$(nproc)
    python test.py
"""
import importlib
import os
import sys
import numpy as np

try:
    import correlation
except ImportError:
    correlation = importlib.import_module("_correlation")


# Force UTF-8 stdout on Windows (cp1252 can't encode box-drawing chars)
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    try:
        reconf_out = getattr(sys.stdout, "reconfigure", None)
        if reconf_out is not None:
            reconf_out(encoding="utf-8", errors="replace")
        reconf_err = getattr(sys.stderr, "reconfigure", None)
        if reconf_err is not None:
            reconf_err(encoding="utf-8", errors="replace")
    except (AttributeError, OSError):
        import io
        sys.stdout = io.TextIOWrapper(
            sys.stdout.buffer, encoding="utf-8", errors="replace"
        )
        sys.stderr = io.TextIOWrapper(
            sys.stderr.buffer, encoding="utf-8", errors="replace"
        )

SEPARATOR = "─" * 60

def section(title):
    print(f"\n{SEPARATOR}")
    print(f"  {title}")
    print(SEPARATOR)

# ── 1. Module import ────────────────────────────────────────────────
section("1. Module")
print(f"  correlation module: OK")
print(f"  KernelType.Gaussian  = {correlation.KernelType.Gaussian}")
print(f"  KernelType.Bump      = {correlation.KernelType.Bump}")
print(f"  KernelType.Triweight = {correlation.KernelType.Triweight}")

# ── 2. Core types ───────────────────────────────────────────────────
section("2. Core types (Atom / Cell / Trajectory)")
cell = correlation.Cell()
a1 = cell.add_atom("Si", [0.0, 0.0, 0.0])
a2 = cell.add_atom("O",  [1.6, 0.0, 0.0])
print(f"  Atoms: {len(cell.atoms)}")
for atom in cell.atoms:
    print(f"    id={atom.id}  element={atom.element.symbol}  pos={atom.position}")

traj = correlation.Trajectory()
print(f"  Trajectory frames: {traj.num_frames()}")

# ── 2.5 NumPy zero-copy bindings ─────────────────────────────────────
section("2.5 NumPy Zero-Copy Bindings")
HAS_NUMPY = False
try:
    import numpy as np
    HAS_NUMPY = True
    
    # Generate some mock data
    pos = cell.positions
    print(f"  cell.positions shape: {pos.shape}, type: {type(pos)}")
    assert pos.shape == (2, 3), "Positions shape mismatch"
    
    # Modify numpy array in-place and verify zero-copy behavior
    pos[0, 0] = 9.5
    assert abs(cell.atoms[0].position[0] - 9.5) < 1e-4, "Zero-copy modification failed!"
    print("  Zero-copy positions verified.")
    
    # Verify velocities
    vel = cell.velocities
    print(f"  cell.velocities shape: {vel.shape}")
    assert vel.shape == (2, 3), "Velocities shape mismatch"
    vel[1, 2] = -4.5
    # We can't access atom.velocity from python if we didn't bind it, but we can read the array again
    assert cell.velocities[1, 2] == -4.5, "Velocity update failed"
    print("  Zero-copy velocities verified.")
except ImportError:
    print("  numpy not installed. Skipping zero-copy test.")

# ── 3. IO (file reading) ─────────────────────────────────────────────
section("3. IO — read()")
try:
    correlation.read("nonexistent.xyz")
    assert False, "Should have thrown RuntimeError for nonexistent file"
except RuntimeError as e:
    print(f"  Expected error (no reader for .xyz or file missing): {e}")

# ── 4. AnalysisSettings ──────────────────────────────────────────────
section("4. AnalysisSettings")
settings = correlation.AnalysisSettings()
settings.r_max       = 15.0
settings.r_bin_width = 0.02
settings.smoothing   = True
settings.smoothing_sigma  = 0.05
settings.smoothing_kernel = correlation.KernelType.Gaussian
settings.active_calculators = {"RDF": True}
print(f"  r_max        = {settings.r_max}")
print(f"  r_bin_width  = {settings.r_bin_width}")
print(f"  smoothing    = {settings.smoothing}")
print(f"  RDF active   = {settings.is_active('RDF')}")
print(f"  SQ active    = {settings.is_active('SQ')} (not in map -> False)")

# ── 5. DistributionFunctions on a minimal cell ───────────────────────
# ── 5. DistributionFunctions on a trailing / single cell ───────────────────
section("5. DistributionFunctions — single cell")
cell2 = correlation.Cell()
cell2.add_atom("Si", [0.0, 0.0, 0.0])
cell2.add_atom("O",  [1.6, 0.0, 0.0])
cell2.add_atom("O",  [0.0, 1.6, 0.0])

# Resolve trajectory path relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))
xdatcar_path = os.path.join(script_dir, "data", "xdatcar", "Si.xdatcar")

traj2 = correlation.read(xdatcar_path)
assert traj2.num_frames() > 0, "No frames loaded"
assert len(traj2) == traj2.num_frames(), "Trajectory __len__ mismatch"

# Test lazy indexing (__getitem__)
cell_xdat = traj2[0]
assert type(cell_xdat) == correlation.Cell, "traj2[0] type mismatch"
assert len(cell_xdat.atoms) > 0, "Cell has no atoms"

# Test negative indexing
cell_last = traj2[-1]
assert type(cell_last) == correlation.Cell, "traj2[-1] type mismatch"

# Test out of bounds indexing raises IndexError
try:
    traj2[len(traj2)]
    assert False, "traj2[len(traj2)] should have raised IndexError"
except IndexError:
    pass

try:
    traj2[-len(traj2) - 1]
    assert False, "traj2[-len(traj2) - 1] should have raised IndexError"
except IndexError:
    pass

df = correlation.DistributionFunctions(cell_xdat, cutoff=5.0, bond_cutoffs=[[3.0]])
df.calculate_rdf(r_max=5.0, bin_width=0.05)
available = df.get_available_histograms()
print(f"  Available histograms: {available}")

if "g_r" in available:
    h = df.get_histogram("g_r")
    print(f"  g_r bins (first 5): {h.bins[:5]}")
    if "Total" in h.partials:
        print(f"  g_r Total (first 5): {h.partials['Total'][:5]}")

df.smooth_all(sigma=0.05)
print("  smooth_all: OK")


# ── 5.5 Dynamic properties getters ──────────────────────────────────
section("5.5 Dynamic Properties Getters")
df_props = correlation.DistributionFunctions(cell2, 0.0, [])
print(f"  Initial MSD diffusion: {df_props.get_diffusion_coefficient_msd()}")
print(f"  Initial VACF diffusion: {df_props.get_diffusion_coefficient_vacf()}")
print(f"  Initial relaxation time: {df_props.get_relaxation_time()}")
print(f"  Initial Deborah number: {df_props.get_deborah_number()}")
assert df_props.get_diffusion_coefficient_msd() == 0.0
assert df_props.get_diffusion_coefficient_vacf() == 0.0
assert df_props.get_relaxation_time() == 0.0
assert df_props.get_deborah_number() == 0.0
print("  Dynamic properties bindings verified: OK")

# ── 5.6 Non-physical parameter guards in Python ──────────────────────
section("5.6 Non-physical Parameter Guards in Python")
try:
    df_props.calculate_rdf(r_max=-1.0, bin_width=0.05)
    assert False, "Should have thrown ValueError/RuntimeError for negative r_max"
except (ValueError, RuntimeError) as e:
    print(f"  Passed negative r_max test (threw expected error): {e}")

try:
    df_props.calculate_pad(bin_width=-1.0)
    assert False, "Should have thrown ValueError/RuntimeError for negative bin_width"
except (ValueError, RuntimeError) as e:
    print(f"  Passed negative bin_width test (threw expected error): {e}")

# ── 6. Calculator access ─────────────────────────────────────────────
section("6. Calculators")
calcs = correlation.get_all_calculators()
print(f"  Registered calculators ({len(calcs)}):")
for c in calcs:
    print(f"    [{c.get_group():10s}] {c.get_short_name():8s} -- {c.get_name()}")

names = correlation.list_calculators()
print(f"  Short names: {names}")

# ── 7.5. GNN Periodic Graph & MLIP ────────────────────────────────────
section("7.5. GNN Periodic Graph & MLIP")
cell_sc = correlation.Cell([4.0, 4.0, 4.0, 90.0, 90.0, 90.0])
cell_sc.add_atom("Si", [0.0, 0.0, 0.0])
# Test PeriodicGraphBuilder.build_graph (cutoff 4.1 encompasses 6 periodic neighbors)
graph = correlation.build_periodic_graph(cell_sc, cutoff_radius=4.1, include_self_loops=False)
print(f"  Graph atoms: {graph.atom_count}, edges: {graph.edge_count}")
assert graph.atom_count == 1, "Graph atom_count mismatch"
assert graph.edge_count == 6, f"Graph edge_count expected 6, got {graph.edge_count}"
if HAS_NUMPY:
    print(f"  Graph pos shape: {graph.positions.shape}")
    print(f"  Graph Z shape: {graph.atomic_numbers.shape}, Z: {graph.atomic_numbers}")
    print(f"  Graph edge_index shape: {graph.edge_index.shape}")
    print(f"  Graph edge_shifts shape: {graph.edge_shifts.shape}")
    print(f"  Graph edge_vectors shape: {graph.edge_vectors.shape}")
    print(f"  Graph edge_distances shape: {graph.edge_distances.shape}")
    assert graph.edge_index.shape == (2, 6), "Edge index shape mismatch"
    assert graph.edge_vectors.shape == (6, 3), "Edge vectors shape mismatch"
    assert graph.edge_distances.shape == (6,), "Edge distances shape mismatch"
    assert all(abs(d - 4.0) < 1e-4 for d in graph.edge_distances), "Edge distance value mismatch"
else:
    print("  numpy not installed. Skipping zero-copy graph array shape checks.")

# Test RBF utilities
env = correlation.PeriodicGraphBuilder.compute_cutoff_envelope(2.5, 5.0)
print(f"  Cutoff envelope(2.5, 5.0): {env}")
assert 0.0 < env < 1.0, "Cutoff envelope out of range"

bessel = correlation.PeriodicGraphBuilder.compute_bessel_basis(2.5, 5.0, 6)
print(f"  Bessel basis (6): {bessel}")
assert len(bessel) == 6, "Bessel basis length mismatch"

rbf = correlation.PeriodicGraphBuilder.compute_gaussian_rbf(1.0, 0.0, 4.0, 5)
print(f"  Gaussian RBF (5): {rbf}")
assert len(rbf) == 5, "Gaussian RBF length mismatch"
assert abs(rbf[1] - 1.0) < 1e-4, "Gaussian RBF peak center mismatch"

# Test MLIPCalculator
mlip_calc = correlation.MLIPCalculator()
print(f"  MLIPCalculator: {mlip_calc.get_name()}")
mlip_out = correlation.MLIPCalculator.calculate(cell_sc)
if HAS_NUMPY:
    print(f"  MLIP output total_energy: {mlip_out.total_energy}, forces shape: {mlip_out.forces.shape}")
    assert mlip_out.forces.shape == (1, 3), "MLIP forces shape mismatch"
else:
    print(f"  MLIP output total_energy: {mlip_out.total_energy}")

# Test TDOSCalculator & TDOSParams
tdos_calc = correlation.TDOSCalculator()
print(f"  TDOSCalculator: {tdos_calc.get_name()} ({tdos_calc.get_short_name()})")
tdos_params = correlation.TDOSParams()
assert tdos_params.e_min == -15.0, "TDOSParams e_min mismatch"
assert tdos_params.e_max == 5.0, "TDOSParams e_max mismatch"
tdos_hist = correlation.TDOSCalculator.calculate(cell_sc, tdos_params)
assert len(tdos_hist.bins) == 0, "TDOS without model should return empty bins"
print("  TDOSCalculator & TDOSParams validated OK")

# ── 7.6. ORB-v3 Graph Descriptors ──────────────────────────────────
section("7.6. ORB-v3 Native Graph Descriptors")
orb_cfg = correlation.OrbDescriptorConfig(r_max=6.0, num_rbf=8, l_max=3, include_self_loops=False, compute_orb_features=True)
assert orb_cfg.r_max == 6.0, "OrbDescriptorConfig r_max mismatch"
assert orb_cfg.num_rbf == 8, "OrbDescriptorConfig num_rbf mismatch"
assert orb_cfg.l_max == 3, "OrbDescriptorConfig l_max mismatch"

orb_graph = correlation.build_orb_graph(cell_sc, orb_cfg)
print(f"  ORB Graph atoms: {orb_graph.atom_count}, edges: {orb_graph.edge_count}")
assert orb_graph.atom_count == 1, "ORB graph atom count mismatch"
assert orb_graph.edge_count > 0, "ORB graph edge count should be > 0"

if HAS_NUMPY:
    print(f"  ORB edge_unit_vectors shape: {orb_graph.edge_unit_vectors.shape}")
    print(f"  ORB edge_radial_basis shape: {orb_graph.edge_radial_basis.shape}")
    print(f"  ORB edge_spherical_harmonics shape: {orb_graph.edge_spherical_harmonics.shape}")
    print(f"  ORB edge_cutoff_envelope shape: {orb_graph.edge_cutoff_envelope.shape}")
    print(f"  ORB edge_orb_features shape: {orb_graph.edge_orb_features.shape}")

    E = orb_graph.edge_count
    assert orb_graph.edge_unit_vectors.shape == (E, 3), "edge_unit_vectors shape mismatch"
    assert orb_graph.edge_radial_basis.shape == (E, 8), "edge_radial_basis shape mismatch"
    assert orb_graph.edge_spherical_harmonics.shape == (E, 16), "edge_spherical_harmonics shape mismatch"
    assert orb_graph.edge_cutoff_envelope.shape == (E,), "edge_cutoff_envelope shape mismatch"
    assert orb_graph.edge_orb_features.shape == (E, 128), "edge_orb_features shape mismatch"

    # Verify unit vector normalization
    norms = np.linalg.norm(orb_graph.edge_unit_vectors, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-4), "Unit vectors not normalized"

    # Verify fused feature values: cutoff * RBF * SH
    c = orb_graph.edge_cutoff_envelope[0]
    r = orb_graph.edge_radial_basis[0]
    y = orb_graph.edge_spherical_harmonics[0]
    fused = orb_graph.edge_orb_features[0]
    expected_fused = (c * np.outer(r, y)).reshape(-1)
    assert np.allclose(fused, expected_fused, atol=1e-4), "Fused ORB features mismatch"
    print("  ORB fused features verification passed.")

# Test standalone math functions
orb_cut = correlation.PeriodicGraphBuilder.compute_orb_cutoff_envelope(3.0, 6.0)
assert abs(orb_cut - 0.65625) < 1e-5, f"Cutoff midpoint mismatch: {orb_cut}"
print(f"  ORB cutoff envelope midpoint: {orb_cut}")

orb_bessel = correlation.PeriodicGraphBuilder.compute_orb_bessel_basis(3.0, orb_cfg)
assert len(orb_bessel) == 8, "Bessel basis length mismatch"
print(f"  ORB Bessel basis: {orb_bessel}")

orb_sh = correlation.PeriodicGraphBuilder.compute_orb_spherical_harmonics([0.0, 1.0, 0.0], 3)
assert len(orb_sh) == 16, "Spherical harmonics length mismatch"
if HAS_NUMPY:
    sh_arr = np.array(orb_sh)
    # On y-axis, sum of squares per l degree should equal 2l+1
    assert abs(sh_arr[0]**2 - 1.0) < 1e-4, "l=0 norm mismatch"
    assert abs(np.sum(sh_arr[1:4]**2) - 3.0) < 1e-4, "l=1 norm mismatch"
    assert abs(np.sum(sh_arr[4:9]**2) - 5.0) < 1e-4, "l=2 norm mismatch"
    assert abs(np.sum(sh_arr[9:16]**2) - 7.0) < 1e-4, "l=3 norm mismatch"
print("  ORB standalone math functions verified OK")

# ── 7.7. Topological Graph Descriptors & Structural-Electronic Correlation ─────
section("7.7. Topological Descriptors & Structural-Electronic Correlation")

fcc_cell = correlation.Cell([3.615, 3.615, 3.615, 90.0, 90.0, 90.0])
fcc_cell.add_atom("Cu", [0.0, 0.0, 0.0])
fcc_cell.add_atom("Cu", [0.0, 1.8075, 1.8075])
fcc_cell.add_atom("Cu", [1.8075, 0.0, 1.8075])
fcc_cell.add_atom("Cu", [1.8075, 1.8075, 0.0])

fcc_graph = correlation.build_periodic_graph(fcc_cell, 2.8)
coords = correlation.compute_coordination_embedding(fcc_graph)
assert len(coords) == 4, "Coordination array length mismatch"
print(f"  FCC Coordination: {coords}")

rings = correlation.compute_ring_statistics_descriptor(fcc_graph, 4)
assert len(rings) == 4 * 4, "Ring statistics descriptor length mismatch"
print(f"  FCC Ring descriptors (N x 4): {len(rings)} elements")

spectrum = correlation.compute_graph_spectrum(fcc_graph, 3)
assert len(spectrum) == 3, "Graph spectrum length mismatch"
assert spectrum[0] >= spectrum[1] >= spectrum[2], "Spectrum not sorted descending"
print(f"  Graph Spectrum top-3 eigenvalues: {spectrum}")

correlation.populate_descriptors(fcc_graph, 4)
if HAS_NUMPY:
    assert len(fcc_graph.cna_labels) == 4
    assert len(fcc_graph.coordination_desc) == 4
    assert fcc_graph.ring_desc.shape == (4, 4)
    assert fcc_graph.ring_desc.size == 16
print("  populate_descriptors verified OK")

# StructuralElectronicCorrelation verification
fcc_traj = correlation.Trajectory()
fcc_traj.add_frame(fcc_cell)
tdos_dists = correlation.DistributionFunctions(fcc_cell)
tdos_params = correlation.TDOSParams(-10.0, 5.0, None)

motif_cna = correlation.correlate_cna(tdos_dists, fcc_traj, tdos_params)
assert isinstance(motif_cna, correlation.MotifProjectedTDOS)
motif_steinhardt = correlation.correlate_steinhardt(tdos_dists, fcc_traj, tdos_params)
assert isinstance(motif_steinhardt, correlation.MotifProjectedTDOS)
print("  StructuralElectronicCorrelation bindings verified OK")

# ── 8. Summary ───────────────────────────────────────────────────────
section("Summary")
print("  All binding layers loaded successfully OK")
print(f"  Python {sys.version}")
