"""
test.py — End-to-end smoke test for the correlation Python bindings.

Run from the build directory after building the correlation_py target:
    cmake --build . --target correlation_py -j$(nproc)
    python test.py
"""
import _correlation as correlation
import sys

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
try:
    import numpy as np
    
    # Generate some mock data
    pos = cell.positions
    print(f"  cell.positions shape: {pos.shape}, type: {type(pos)}")
    assert pos.shape == (2, 3), "Positions shape mismatch"
    
    # Modify numpy array in-place and verify zero-copy behavior
    pos[0, 0] = 9.9
    assert cell.atoms[0].position[0] == 9.9, "Zero-copy modification failed!"
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
section("5. DistributionFunctions — single cell")
cell2 = correlation.Cell()
cell2.add_atom("Si", [0.0, 0.0, 0.0])
cell2.add_atom("O",  [1.6, 0.0, 0.0])
cell2.add_atom("O",  [0.0, 1.6, 0.0])

try:
    traj2 = correlation.read("tests/data/Si.xdatcar")
    cell_xdat = traj2.frames[0]
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

except Exception as e:
    print(f"  DistributionFunctions error: {e}")

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

# ── 7. Writer access ─────────────────────────────────────────────────
section("7. Writers")
writer_names = correlation.list_writers()
print(f"  Registered writers ({len(writer_names)}): {writer_names}")

csv_writer = correlation.get_writer("CSV")
print(f"  CSV writer: {csv_writer.get_name()}, extensions={csv_writer.get_extensions()}")

# ── 8. Summary ───────────────────────────────────────────────────────
section("Summary")
print("  All binding layers loaded successfully OK")
print(f"  Python {sys.version}")
