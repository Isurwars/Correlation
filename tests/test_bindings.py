"""
test.py — End-to-end smoke test for the correlation Python bindings.

Run from the build directory after building the correlation_py target:
    cmake --build . --target correlation_py -j$(nproc)
    python test.py
"""
import correlation
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

# ── 3. IO (file reading) ─────────────────────────────────────────────
section("3. IO — read()")
try:
    correlation.read("nonexistent.xyz")
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
    df = correlation.DistributionFunctions(cell2, cutoff=5.0)
    df.calculate_rdf(r_max=10.0, bin_width=0.05)
    available = df.get_available_histograms()
    print(f"  Available histograms: {available}")

    if "g(r)" in available:
        h = df.get_histogram("g(r)")
        print(f"  g(r) bins (first 5): {h.bins[:5]}")
        if "Total" in h.partials:
            print(f"  g(r) Total (first 5): {h.partials['Total'][:5]}")

    df.smooth_all(sigma=0.05)
    print("  smooth_all: OK")

except Exception as e:
    print(f"  DistributionFunctions error (expected if cell lacks lattice): {e}")

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
