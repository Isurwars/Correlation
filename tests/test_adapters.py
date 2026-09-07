"""
Unit tests for correlation.adapters: ASE and Pymatgen ecosystem bridges.
"""

import sys
from types import ModuleType
from unittest.mock import MagicMock
import numpy as np
import pytest

import correlation
from correlation import adapters


# -----------------------------------------------------------------------------
# Fixtures for Mocking ASE and Pymatgen
# -----------------------------------------------------------------------------

@pytest.fixture
def mock_ase():
    """Mock the ase module and ase.Atoms class."""
    mock_mod = ModuleType("ase")

    class MockCell:
        def __init__(self, cellpar):
            self._cellpar = list(cellpar) if cellpar is not None else [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            self.volume = 1000.0 if any(p > 0 for p in self._cellpar) else 0.0

        def cellpar(self):
            return self._cellpar

    class MockAtoms:
        def __init__(self, symbols=None, positions=None, cell=None, pbc=False):
            self.symbols = list(symbols) if symbols is not None else []
            self.positions = np.array(positions, dtype=float) if positions is not None else np.empty((0, 3))
            self.pbc = [pbc, pbc, pbc] if isinstance(pbc, bool) else list(pbc)
            self.cell = MockCell(cell) if cell is not None else MockCell([0, 0, 0, 0, 0, 0])

        def get_chemical_symbols(self):
            return self.symbols

        def get_positions(self):
            return self.positions

    mock_mod.Atoms = MockAtoms
    sys.modules["ase"] = mock_mod
    yield mock_mod
    sys.modules.pop("ase", None)


@pytest.fixture
def mock_pymatgen():
    """Mock pymatgen and pymatgen.core modules."""
    pmg_mod = ModuleType("pymatgen")
    pmg_core = ModuleType("pymatgen.core")

    class MockLattice:
        def __init__(self, a, b, c, alpha, beta, gamma):
            self.a = float(a)
            self.b = float(b)
            self.c = float(c)
            self.alpha = float(alpha)
            self.beta = float(beta)
            self.gamma = float(gamma)

        @classmethod
        def from_parameters(cls, a, b, c, alpha, beta, gamma):
            return cls(a, b, c, alpha, beta, gamma)

    class MockSpecie:
        def __init__(self, symbol):
            self.symbol = symbol

        def __str__(self):
            return self.symbol

    class MockSite:
        def __init__(self, specie, coords):
            self.specie = MockSpecie(specie) if isinstance(specie, str) else specie
            self.coords = np.array(coords, dtype=float)

    class MockStructure:
        def __init__(self, lattice, species, coords, coords_are_cartesian=True):
            self.lattice = lattice
            self.sites = [MockSite(s, c) for s, c in zip(species, coords)]

        def __iter__(self):
            return iter(self.sites)

        def __len__(self):
            return len(self.sites)

    class MockMolecule:
        def __init__(self, species, coords):
            self.sites = [MockSite(s, c) for s, c in zip(species, coords)]

        def __iter__(self):
            return iter(self.sites)

        def __len__(self):
            return len(self.sites)

    pmg_core.Lattice = MockLattice
    pmg_core.Structure = MockStructure
    pmg_core.Molecule = MockMolecule
    pmg_core.Site = MockSite
    pmg_mod.core = pmg_core

    sys.modules["pymatgen"] = pmg_mod
    sys.modules["pymatgen.core"] = pmg_core
    yield pmg_core
    sys.modules.pop("pymatgen.core", None)
    sys.modules.pop("pymatgen", None)


# -----------------------------------------------------------------------------
# Tests: Dependency Absence
# -----------------------------------------------------------------------------

def test_ase_missing_raises_import_error():
    """Verify informative ImportError when ASE is not installed."""
    sys.modules.pop("ase", None)
    cell = correlation.Cell([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0])

    with pytest.raises(ImportError, match="Atomic Simulation Environment"):
        adapters.from_ase(None)

    with pytest.raises(ImportError, match="Atomic Simulation Environment"):
        adapters.to_ase(cell)

    with pytest.raises(ImportError, match="Atomic Simulation Environment"):
        adapters.from_ase_trajectory([])

    traj = correlation.Trajectory()
    traj.add_frame(cell)
    with pytest.raises(ImportError, match="Atomic Simulation Environment"):
        adapters.to_ase_trajectory(traj)


def test_pymatgen_missing_raises_import_error():
    """Verify informative ImportError when Pymatgen is not installed."""
    sys.modules.pop("pymatgen", None)
    sys.modules.pop("pymatgen.core", None)
    cell = correlation.Cell([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0])

    with pytest.raises(ImportError, match="Pymatgen is required"):
        adapters.from_pymatgen(None)

    with pytest.raises(ImportError, match="Pymatgen is required"):
        adapters.to_pymatgen(cell)


# -----------------------------------------------------------------------------
# Tests: ASE Conversion Logic
# -----------------------------------------------------------------------------

def test_from_ase_periodic(mock_ase):
    """Test converting periodic ASE Atoms to Correlation Cell."""
    atoms = mock_ase.Atoms(
        symbols=["Si", "O", "O"],
        positions=[[0.0, 0.0, 0.0], [1.2, 1.2, 0.0], [2.4, 0.0, 1.2]],
        cell=[10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
        pbc=True,
    )
    cell = adapters.from_ase(atoms)

    assert len(cell) == 3
    assert [a.element.symbol for a in cell.atoms] == ["Si", "O", "O"]
    np.testing.assert_allclose(cell.positions, [[0.0, 0.0, 0.0], [1.2, 1.2, 0.0], [2.4, 0.0, 1.2]], rtol=1e-5)
    params = cell.get_lattice_parameters()
    np.testing.assert_allclose(params[:3], [10.0, 10.0, 10.0], rtol=1e-4)


def test_from_ase_non_periodic(mock_ase):
    """Test converting isolated/cluster ASE Atoms to non-periodic Correlation Cell."""
    atoms = mock_ase.Atoms(
        symbols=["H", "H", "O"],
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]],
        cell=None,
        pbc=False,
    )
    cell = adapters.from_ase(atoms)

    assert len(cell) == 3
    assert [a.element.symbol for a in cell.atoms] == ["H", "H", "O"]
    np.testing.assert_allclose(cell.positions, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]], rtol=1e-5)


def test_to_ase_periodic(mock_ase):
    """Test converting periodic Correlation Cell to ASE Atoms."""
    cell = correlation.Cell([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0])
    cell.add_atom("C", [1.0, 1.0, 1.0])
    cell.add_atom("O", [2.0, 2.0, 2.0])

    atoms = adapters.to_ase(cell)
    assert atoms.get_chemical_symbols() == ["C", "O"]
    np.testing.assert_allclose(atoms.get_positions(), [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]], rtol=1e-5)
    assert all(atoms.pbc)
    np.testing.assert_allclose(atoms.cell.cellpar()[:3], [10.0, 10.0, 10.0], rtol=1e-4)


def test_to_ase_non_periodic(mock_ase):
    """Test converting empty/cluster Cell to ASE Atoms with pbc=False."""
    cell = correlation.Cell()
    cell.add_atom("He", [0.0, 0.0, 0.0])

    atoms = adapters.to_ase(cell)
    assert atoms.get_chemical_symbols() == ["He"]
    assert not any(atoms.pbc)


def test_ase_trajectory_roundtrip(mock_ase):
    """Test multi-frame trajectory conversion between ASE and Correlation."""
    atoms1 = mock_ase.Atoms(symbols=["Ar"], positions=[[0.0, 0.0, 0.0]], cell=[5.0, 5.0, 5.0, 90, 90, 90], pbc=True)
    atoms2 = mock_ase.Atoms(symbols=["Ar"], positions=[[0.5, 0.5, 0.5]], cell=[5.0, 5.0, 5.0, 90, 90, 90], pbc=True)

    traj = adapters.from_ase_trajectory([atoms1, atoms2])
    assert len(traj) == 2
    assert len(traj[0]) == 1
    assert len(traj[1]) == 1

    converted_images = adapters.to_ase_trajectory(traj)
    assert len(converted_images) == 2
    assert converted_images[0].get_chemical_symbols() == ["Ar"]
    assert converted_images[1].get_chemical_symbols() == ["Ar"]


# -----------------------------------------------------------------------------
# Tests: Pymatgen Conversion Logic
# -----------------------------------------------------------------------------

def test_from_pymatgen_structure(mock_pymatgen):
    """Test converting Pymatgen Structure to Correlation Cell."""
    lattice = mock_pymatgen.Lattice.from_parameters(8.0, 8.0, 8.0, 90.0, 90.0, 90.0)
    structure = mock_pymatgen.Structure(lattice, ["Fe", "Fe"], [[0.0, 0.0, 0.0], [4.0, 4.0, 4.0]])

    cell = adapters.from_pymatgen(structure)
    assert len(cell) == 2
    assert [a.element.symbol for a in cell.atoms] == ["Fe", "Fe"]
    np.testing.assert_allclose(cell.positions, [[0.0, 0.0, 0.0], [4.0, 4.0, 4.0]], rtol=1e-5)


def test_from_pymatgen_molecule(mock_pymatgen):
    """Test converting Pymatgen Molecule to Correlation Cell."""
    molecule = mock_pymatgen.Molecule(["H", "F"], [[0.0, 0.0, 0.0], [0.92, 0.0, 0.0]])

    cell = adapters.from_pymatgen(molecule)
    assert len(cell) == 2
    assert [a.element.symbol for a in cell.atoms] == ["H", "F"]
    np.testing.assert_allclose(cell.positions, [[0.0, 0.0, 0.0], [0.92, 0.0, 0.0]], rtol=1e-5)


def test_to_pymatgen_structure(mock_pymatgen):
    """Test converting periodic Cell to Pymatgen Structure."""
    cell = correlation.Cell([7.0, 0.0, 0.0], [0.0, 7.0, 0.0], [0.0, 0.0, 7.0])
    cell.add_atom("Na", [0.0, 0.0, 0.0])
    cell.add_atom("Cl", [3.5, 3.5, 3.5])

    structure = adapters.to_pymatgen(cell)
    assert isinstance(structure, mock_pymatgen.Structure)
    assert len(structure) == 2
    assert [site.specie.symbol for site in structure] == ["Na", "Cl"]


def test_to_pymatgen_molecule(mock_pymatgen):
    """Test converting non-periodic Cell to Pymatgen Molecule."""
    cell = correlation.Cell()
    cell.add_atom("N", [0.0, 0.0, 0.0])
    cell.add_atom("N", [1.1, 0.0, 0.0])

    molecule = adapters.to_pymatgen(cell)
    assert isinstance(molecule, mock_pymatgen.Molecule)
    assert len(molecule) == 2
    assert [site.specie.symbol for site in molecule] == ["N", "N"]


# -----------------------------------------------------------------------------
# Tests: Monkey-Patched Methods on Cell and Trajectory
# -----------------------------------------------------------------------------

def test_monkey_patched_methods(mock_ase, mock_pymatgen):
    """Verify to_ase, from_ase, to_pymatgen, from_pymatgen work via Cell and Trajectory."""
    cell = correlation.Cell([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0])
    cell.add_atom("Si", [0.0, 0.0, 0.0])

    # Cell.to_ase()
    atoms = cell.to_ase()
    assert atoms.get_chemical_symbols() == ["Si"]

    # Cell.from_ase()
    cell_from_ase = correlation.Cell.from_ase(atoms)
    assert len(cell_from_ase) == 1
    assert cell_from_ase.atoms[0].element.symbol == "Si"

    # Cell.to_pymatgen()
    structure = cell.to_pymatgen()
    assert len(structure) == 1

    # Cell.from_pymatgen()
    cell_from_pmg = correlation.Cell.from_pymatgen(structure)
    assert len(cell_from_pmg) == 1

    # Trajectory.to_ase() and Trajectory.from_ase()
    traj = correlation.Trajectory()
    traj.add_frame(cell)
    ase_traj = traj.to_ase()
    assert len(ase_traj) == 1

    traj_from_ase = correlation.Trajectory.from_ase(ase_traj)
    assert len(traj_from_ase) == 1
