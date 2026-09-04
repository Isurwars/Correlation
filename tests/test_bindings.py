"""
test.py — End-to-end smoke test for the correlation Python bindings.

Run from the build directory after building the correlation_py target:
    cmake --build . --target correlation_py -j$(nproc)
    python test.py
"""
import importlib
import os
import sys

try:
    import correlation
except ImportError:
    correlation = importlib.import_module("_correlation")


# Force UTF-8 stdout on Windows (cp1252 can't encode box-drawing chars)
if sys.stdout.encoding and sys.stdout.encoding.lower() not in ("utf-8", "utf8"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
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

# ── 8. Summary ───────────────────────────────────────────────────────
section("Summary")
print("  All binding layers loaded successfully OK")
print(f"  Python {sys.version}")
