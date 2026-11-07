"""
Unit tests for ASEMLIPBridge, MaceGNNModel bindings, and ASE ecosystem interoperability.
"""

import sys
from types import ModuleType
import numpy as np
import pytest

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import correlation
    from correlation import ase_bridge, adapters
else:
    try:
        import matcorr as correlation
        from matcorr import ase_bridge, adapters
    except (ImportError, AttributeError):
        try:
            import correlation
            from correlation import ase_bridge, adapters
        except (ImportError, AttributeError):
            correlation = None
            ase_bridge = None
            adapters = None

if correlation is None or not hasattr(correlation, "Cell"):
    pytest.skip(
        "correlation/matcorr C++ extension module is not installed (run 'pip install .' to test bindings)",
        allow_module_level=True,
    )


# -----------------------------------------------------------------------------
# Fixtures for Mocking ASE
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
            self.calc = None

        def __len__(self):
            return len(self.symbols)

        def get_chemical_symbols(self):
            return self.symbols

        def get_positions(self):
            return self.positions

        def get_potential_energy(self):
            if self.calc:
                return self.calc.get_potential_energy(self)
            return 0.0

        def get_forces(self):
            if self.calc:
                return self.calc.get_forces(self)
            return np.zeros((len(self.positions), 3))

        def get_stress(self, voigt=False):
            if self.calc:
                return self.calc.get_stress(self, voigt=voigt)
            return np.zeros((3, 3)) if not voigt else np.zeros(6)

        def get_charges(self):
            if self.calc and hasattr(self.calc, "get_charges"):
                return self.calc.get_charges(self)
            return None

        def get_magnetic_moments(self):
            if self.calc and hasattr(self.calc, "get_magnetic_moments"):
                return self.calc.get_magnetic_moments(self)
            return None

    setattr(mock_mod, "Atoms", MockAtoms)
    sys.modules["ase"] = mock_mod
    yield mock_mod
    sys.modules.pop("ase", None)


# -----------------------------------------------------------------------------
# Mock Calculator Fixtures
# -----------------------------------------------------------------------------

class MockMaceCalculator:
    """Mock ASE Calculator mimicking MACE-MP-0 / MACE-OFF."""
    name = "mace_mp"

    def __init__(self, energy_per_atom: float = -5.43):
        self.energy_per_atom = energy_per_atom

    def get_potential_energy(self, atoms=None):
        n = len(atoms) if atoms is not None else 1
        return self.energy_per_atom * n

    def get_forces(self, atoms=None):
        n = len(atoms) if atoms is not None else 1
        # Deterministic dummy forces
        return np.tile([0.1, -0.2, 0.05], (n, 1))

    def get_stress(self, atoms=None, voigt=False):
        # 3x3 diagonal pressure
        return np.diag([-0.01, -0.01, -0.01])

    def get_charges(self, atoms=None):
        n = len(atoms) if atoms is not None else 1
        return np.ones(n) * 0.5


class MockSevenNetCalculator:
    """Mock ASE Calculator mimicking SevenNet-0."""
    name = "sevennet_0"

    def get_potential_energy(self, atoms=None):
        return -120.5

    def get_forces(self, atoms=None):
        n = len(atoms) if atoms is not None else 1
        return np.zeros((n, 3))

    def get_stress(self, atoms=None, voigt=False):
        # 6-element Voigt stress: [xx, yy, zz, yz, xz, xy]
        return np.array([-0.05, -0.05, -0.05, 0.0, 0.0, 0.0])


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

def test_mace_gnn_model_bindings():
    """Verify MaceGNNModel C++ bindings, accessors, and fallback behavior."""
    if not hasattr(correlation, "MaceGNNModel"):
        pytest.skip("MaceGNNModel is not exposed in bindings.")

    model = correlation.MaceGNNModel("dummy_mace.pt", "cpu", 5.5)
    assert model.cutoff == pytest.approx(5.5)
    model.cutoff = 4.8
    assert model.cutoff == pytest.approx(4.8)
    assert "dummy_mace.pt" in model.get_model_name()
    assert model.get_model_name().startswith("MACE [")

    # Empty cell evaluation
    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    out = model.evaluate(cell)
    assert isinstance(out, correlation.MLIPOutput)
    assert len(out.forces) == 0
    assert out.total_energy == pytest.approx(0.0)


def test_mlip_output_properties():
    """Verify MLIPOutput property getters and setters for forces and stress."""
    out = correlation.MLIPOutput()
    out.total_energy = -42.5
    assert out.total_energy == pytest.approx(-42.5)

    # Set and get forces
    test_forces = np.array([[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]], dtype=np.float64)
    out.forces = test_forces
    forces_ret = out.forces
    assert forces_ret.shape == (2, 3)
    np.testing.assert_allclose(forces_ret, test_forces)

    # Set and get stress
    test_stress = np.eye(3, dtype=np.float64) * 0.123
    out.stress = test_stress
    stress_ret = out.stress
    assert stress_ret.shape == (3, 3)
    np.testing.assert_allclose(stress_ret, test_stress)


def test_ase_mlip_bridge_mace(mock_ase):
    """Test ASEMLIPBridge with mock MACE calculator."""
    calc = MockMaceCalculator()
    bridge = ase_bridge.ASEMLIPBridge(calc)

    assert "mace_mp" in bridge.get_model_name()

    cell = correlation.Cell([12.0, 12.0, 12.0, 90.0, 90.0, 90.0])
    cell.add_atom("Si", [0.0, 0.0, 0.0])
    cell.add_atom("O", [1.6, 0.0, 0.0])

    output = bridge.evaluate(cell)
    assert isinstance(output, correlation.MLIPOutput)
    assert output.total_energy == pytest.approx(-10.86)
    assert len(output.forces) == 2
    np.testing.assert_allclose(output.forces[0], [0.1, -0.2, 0.05])
    assert output.stress.shape == (3, 3)


def test_ase_mlip_bridge_sevennet_voigt(mock_ase):
    """Test ASEMLIPBridge Voigt 6-element stress conversion."""
    calc = MockSevenNetCalculator()
    bridge = ase_bridge.ASEMLIPBridge(calc)

    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    cell.add_atom("Fe", [0.0, 0.0, 0.0])

    output = bridge.evaluate(cell)
    assert output.total_energy == pytest.approx(-120.5)
    assert output.stress.shape == (3, 3)
    assert output.stress[0, 0] == pytest.approx(-0.05)
    assert output.stress[1, 1] == pytest.approx(-0.05)
    assert output.stress[2, 2] == pytest.approx(-0.05)


def test_cell_and_trajectory_evaluate_mlip(mock_ase):
    """Test convenience evaluate_mlip methods attached to Cell and Trajectory."""
    calc = MockMaceCalculator()

    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    cell.add_atom("C", [0.0, 0.0, 0.0])

    out = cell.evaluate_mlip(calc)
    assert isinstance(out, correlation.MLIPOutput)
    assert out.total_energy == pytest.approx(-5.43)

    traj = correlation.Trajectory()
    traj.add_frame(cell)
    traj.add_frame(cell)

    traj_results = traj.evaluate_mlip(calc)
    assert len(traj_results) == 2
    assert traj_results[0].total_energy == pytest.approx(-5.43)
    assert traj_results[1].total_energy == pytest.approx(-5.43)


def test_cpp_mlip_calculator_with_ase_bridge(mock_ase):
    """Test passing Python ASEMLIPBridge directly into C++ MLIPCalculator."""
    calc = MockMaceCalculator()
    bridge = ase_bridge.ASEMLIPBridge(calc)

    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    cell.add_atom("Si", [0.0, 0.0, 0.0])

    if hasattr(correlation, "MLIPCalculator"):
        cxx_out = correlation.MLIPCalculator.calculate(cell, bridge)
        assert isinstance(cxx_out, correlation.MLIPOutput)
        assert cxx_out.total_energy == pytest.approx(-5.43)
        assert len(cxx_out.forces) == 1


def test_spin_correlation_calculator():
    """Verify SpinCorrelationCalculator on ferromagnetically and antiferromagnetically aligned spins."""
    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    cell.add_atom("Fe", [0.0, 0.0, 0.0])
    cell.add_atom("Fe", [2.5, 0.0, 0.0])

    # Ferromagnetic alignment: [0, 0, 2.2] and [0, 0, 2.2]
    ferro_spins = np.array([[0.0, 0.0, 2.2], [0.0, 0.0, 2.2]], dtype=np.float64)
    res_ferro = correlation.SpinCorrelationCalculator.calculate(cell, ferro_spins, r_max=5.0, r_bin_width=0.5)

    assert len(res_ferro.r_bins) == 10
    # Bin around 2.5 A (index 5)
    bin_idx = int(2.5 / 0.5)
    assert res_ferro.pair_counts[bin_idx] == 1
    assert res_ferro.spin_correlation[bin_idx] == pytest.approx(1.0)
    assert res_ferro.magnetic_moment_corr[bin_idx] == pytest.approx(2.2 * 2.2)

    # Antiferromagnetic alignment: [0, 0, 2.2] and [0, 0, -2.2]
    antiferro_spins = np.array([[0.0, 0.0, 2.2], [0.0, 0.0, -2.2]], dtype=np.float64)
    res_anti = correlation.SpinCorrelationCalculator.calculate(cell, antiferro_spins, r_max=5.0, r_bin_width=0.5)
    assert res_anti.spin_correlation[bin_idx] == pytest.approx(-1.0)
    assert res_anti.magnetic_moment_corr[bin_idx] == pytest.approx(-2.2 * 2.2)


def test_charge_correlation_calculator():
    """Verify ChargeCorrelationCalculator g_qq(r) on NaCl-like charge pairs."""
    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    cell.add_atom("Na", [0.0, 0.0, 0.0])
    cell.add_atom("Cl", [2.8, 0.0, 0.0])

    # Opposite charges: +1.0 and -1.0
    charges = np.array([1.0, -1.0], dtype=np.float64)
    res = correlation.ChargeCorrelationCalculator.calculate(cell, charges, r_max=5.0, r_bin_width=0.5)

    assert len(res.r_bins) == 10
    bin_idx = int(2.8 / 0.5)
    assert res.pair_counts[bin_idx] == 1
    # Opposite charges produce negative g_qq(r) reflecting screening attraction
    assert res.g_qq[bin_idx] < 0.0
    assert res.charge_integral[bin_idx] < 0.0


def test_mlip_output_spin_and_charge_integration():
    """Verify calculate_from_mlip on MLIPOutput containing charges and magmoms."""
    cell = correlation.Cell([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
    cell.add_atom("Ni", [0.0, 0.0, 0.0])
    cell.add_atom("O", [2.0, 0.0, 0.0])

    out = correlation.MLIPOutput()
    out.charges = np.array([1.8, -1.8], dtype=np.float64)
    out.magmoms = np.array([[0.0, 0.0, 1.7], [0.0, 0.0, 0.0]], dtype=np.float64)

    spin_res = correlation.SpinCorrelationCalculator.calculate_from_mlip(cell, out, r_max=4.0, r_bin_width=0.5)
    charge_res = correlation.ChargeCorrelationCalculator.calculate_from_mlip(cell, out, r_max=4.0, r_bin_width=0.5)

    assert len(spin_res.r_bins) == 8
    assert len(charge_res.r_bins) == 8
    bin_idx = int(2.0 / 0.5)
    assert charge_res.pair_counts[bin_idx] == 1
    assert charge_res.g_qq[bin_idx] < 0.0
