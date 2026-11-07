"""
ase_bridge.py — Bidirectional bridge between Correlation and ASE Machine Learning Potentials.

Provides ASEMLIPBridge, enabling any ASE Calculator (MACE, SevenNet, CHGNet, MatterSim,
ORB, AIMNet2, etc.) to act as a correlation.MLIPInterface for C++ and Python pipelines,
as well as CorrelationASECalculator to expose Correlation analyses within ASE.
"""

from __future__ import annotations

from typing import Any, Callable, List, Optional
import numpy as np

try:
    import correlation
except ImportError:  # pragma: no cover
    try:
        from correlation import _correlation as correlation  # pyright: ignore[reportMissingModuleSource]
    except ImportError:
        correlation = None  # type: ignore[assignment]

from correlation.adapters import _require_ase, to_ase, from_ase


# -----------------------------------------------------------------------------
# Base Class Selection
# -----------------------------------------------------------------------------

# If correlation.MLIPInterface is available, subclass it for C++ polymorphism;
# otherwise fall back to object.
_BaseInterface = getattr(correlation, "MLIPInterface", object) if correlation is not None else object


# -----------------------------------------------------------------------------
# ASEMLIPBridge: Wrap any ASE Calculator as an MLIPInterface
# -----------------------------------------------------------------------------

class ASEMLIPBridge(_BaseInterface):  # type: ignore[misc]
    """
    Adapter wrapping an arbitrary ASE Calculator into a Correlation MLIPInterface.

    Enables state-of-the-art machine learning potentials (e.g. MACE, SevenNet,
    CHGNet, MatterSim, TorchMD-Net, AIMNet2) to evaluate energy, atomic forces,
    and stress tensors directly over Correlation Cell and Trajectory structures.

    Parameters
    ----------
    calculator : ase.calculators.calculator.Calculator or Any
        Configured ASE Calculator instance.
    name : str, optional
        Custom model descriptor name. If omitted, derived from calculator class.
    """

    def __init__(self, calculator: Any, name: Optional[str] = None) -> None:
        if _BaseInterface is not object:
            super().__init__()
        self._calculator = calculator
        calc_name = getattr(calculator, "name", calculator.__class__.__name__)
        self._name = name or f"ASEBridge [{calc_name}]"

    @property
    def calculator(self) -> Any:
        """Return the underlying ASE Calculator instance."""
        return self._calculator

    def get_model_name(self) -> str:
        """
        Return the model descriptor name.

        Returns
        -------
        str
            Descriptive string identifying the model and bridge.
        """
        return self._name

    def evaluate(self, cell: Any) -> Any:
        """
        Evaluate potential energy, atomic forces, and virial stress on a Cell.

        Parameters
        ----------
        cell : correlation.Cell
            Atomic configuration cell to evaluate.

        Returns
        -------
        correlation.MLIPOutput
            Standardized container with evaluated total_energy, forces, and stress.
        """
        if correlation is None or not hasattr(correlation, "MLIPOutput"):
            raise RuntimeError("Correlation C++ extension module is not loaded.")

        output = correlation.MLIPOutput()

        num_atoms = len(cell.atoms) if hasattr(cell, "atoms") else (len(cell.positions) if hasattr(cell, "positions") and cell.positions is not None else 0)
        if num_atoms == 0:
            return output

        # Convert to ASE Atoms
        atoms = to_ase(cell)
        atoms.calc = self._calculator

        # 1. Total Potential Energy
        try:
            output.total_energy = float(atoms.get_potential_energy())
        except Exception:
            output.total_energy = 0.0

        # 2. Atomic Forces (N, 3)
        try:
            forces = np.ascontiguousarray(atoms.get_forces(), dtype=np.float64)
            if forces.ndim == 2 and forces.shape == (num_atoms, 3):
                output.forces = forces
        except Exception:
            pass

        # 3. Stress Tensor (3, 3)
        if any(getattr(atoms, "pbc", [])):
            try:
                # Request 3x3 tensor (voigt=False)
                stress = atoms.get_stress(voigt=False)
                stress_arr = np.ascontiguousarray(stress, dtype=np.float64)
                if stress_arr.shape == (3, 3):
                    output.stress = stress_arr
                elif stress_arr.size == 6:
                    # Convert 6-element Voigt [xx, yy, zz, yz, xz, xy] to 3x3 matrix
                    s_mat = np.zeros((3, 3), dtype=np.float64)
                    s_mat[0, 0] = stress_arr[0]
                    s_mat[1, 1] = stress_arr[1]
                    s_mat[2, 2] = stress_arr[2]
                    s_mat[1, 2] = s_mat[2, 1] = stress_arr[3]
                    s_mat[0, 2] = s_mat[2, 0] = stress_arr[4]
                    s_mat[0, 1] = s_mat[1, 0] = stress_arr[5]
                    output.stress = s_mat
            except Exception:
                pass

        # 4. Optional Per-Atom Quantities (charges, magnetic moments)
        try:
            charges = atoms.get_charges()
            if charges is not None and len(charges) == num_atoms:
                output.charges = np.ascontiguousarray(charges, dtype=np.float64)
        except Exception:
            pass

        try:
            magmoms = atoms.get_magnetic_moments()
            if magmoms is not None and len(magmoms) == num_atoms:
                mag_arr = np.ascontiguousarray(magmoms, dtype=np.float64)
                if mag_arr.ndim == 1:
                    mag_vecs = np.zeros((num_atoms, 3), dtype=np.float64)
                    mag_vecs[:, 2] = mag_arr
                    output.magmoms = mag_vecs
                elif mag_arr.ndim == 2 and mag_arr.shape == (num_atoms, 3):
                    output.magmoms = mag_arr
        except Exception:
            pass

        return output

    def evaluate_trajectory(
        self,
        traj: Any,
        on_progress: Optional[Callable[[float, str], None]] = None,
    ) -> List[Any]:
        """
        Evaluate a sequence of frames in a Correlation Trajectory.

        Parameters
        ----------
        traj : correlation.Trajectory
            Multi-frame trajectory.
        on_progress : Callable[[float, str], None], optional
            Callback receiving (progress_fraction, status_text).

        Returns
        -------
        list[correlation.MLIPOutput]
            List of evaluated results for each frame.
        """
        total_frames = len(traj)
        results: List[Any] = []
        if total_frames == 0:
            return results

        if on_progress:
            on_progress(0.0, f"Starting MLIP evaluation with {self._name}...")

        for idx in range(total_frames):
            frame = traj[idx]
            res = self.evaluate(frame)
            results.append(res)

            if on_progress and (idx % 10 == 0 or idx == total_frames - 1):
                frac = float(idx + 1) / float(total_frames)
                on_progress(frac, f"Evaluated frame {idx + 1}/{total_frames}")

        return results


# -----------------------------------------------------------------------------
# CorrelationASECalculator: Expose Correlation Analyses into ASE Workflows
# -----------------------------------------------------------------------------

def create_correlation_ase_calculator(
    descriptors: Optional[List[str]] = None,
    rdf_cutoff: float = 10.0,
    rdf_bins: int = 200,
) -> Any:
    """
    Factory creating an ASE Calculator subclass backed by Correlation algorithms.

    Parameters
    ----------
    descriptors : list[str], optional
        List of structural properties to evaluate ("rdf", "cna", "angle").
    rdf_cutoff : float, optional
        Distance cutoff for radial distribution functions.
    rdf_bins : int, optional
        Number of bins for histogram discretization.

    Returns
    -------
    ase.calculators.calculator.Calculator
        Configured ASE Calculator.
    """
    _require_ase()
    from ase.calculators.calculator import Calculator, all_changes

    active_descriptors = descriptors or ["rdf"]

    class CorrelationASECalculator(Calculator):
        implemented_properties = ["energy", "forces"] + active_descriptors

        def __init__(self, **kwargs: Any) -> None:
            super().__init__(**kwargs)
            self.results: dict[str, Any] = {}

        def calculate(
            self,
            atoms: Any = None,
            properties: Optional[List[str]] = None,
            system_changes: Any = all_changes,
        ) -> None:
            super().calculate(atoms, properties, system_changes)
            if self.atoms is None or correlation is None:
                return

            cell = from_ase(self.atoms)

            # Basic placeholder energy/forces for pure structural calculators
            self.results["energy"] = 0.0
            self.results["forces"] = np.zeros((len(self.atoms), 3), dtype=float)

            # Compute requested structural correlations
            if "rdf" in active_descriptors:
                calc = correlation.calculators.RDFCalculator()
                self.results["rdf"] = calc.calculate(cell) if hasattr(calc, "calculate") else None

    return CorrelationASECalculator()

