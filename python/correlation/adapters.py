"""
adapters.py — High-performance ecosystem bridges between Correlation and ASE / Pymatgen.

Provides zero-copy and structured conversion utilities between Correlation's C++
Cell and Trajectory types and atomic simulation data structures in ASE (Atoms, Trajectory)
and Pymatgen (Structure, Molecule).
"""

from __future__ import annotations

from typing import Any, Iterable, List
import numpy as np

try:
    import correlation
except ImportError:  # pragma: no cover
    try:
        from correlation import _correlation as correlation
    except ImportError:
        correlation = None  # type: ignore[assignment]


# -----------------------------------------------------------------------------
# Dependency Check Helpers
# -----------------------------------------------------------------------------

def _require_ase():
    """Verify ASE installation with informative diagnostic."""
    try:
        import ase  # noqa: F401
        from ase import Atoms
        return Atoms
    except ImportError as err:
        raise ImportError(
            "ASE (Atomic Simulation Environment) is required for this operation. "
            "Install it via 'pip install ase'."
        ) from err


def _require_pymatgen():
    """Verify Pymatgen installation with informative diagnostic."""
    try:
        import pymatgen.core as pmg_core  # noqa: F401
        from pymatgen.core import Lattice, Molecule, Structure
        return Structure, Molecule, Lattice
    except ImportError as err:
        raise ImportError(
            "Pymatgen is required for this operation. "
            "Install it via 'pip install pymatgen'."
        ) from err


# -----------------------------------------------------------------------------
# ASE Adapters
# -----------------------------------------------------------------------------

def from_ase(atoms: Any) -> Any:
    """
    Convert an ASE Atoms object into a Correlation Cell.

    Parameters
    ----------
    atoms : ase.Atoms
        Source ASE atomic configuration.

    Returns
    -------
    correlation.Cell
        Constructed Simulation Cell with matching lattice vectors and atomic positions.

    Raises
    ------
    ImportError
        If ASE is not installed in the current environment.
    """
    _require_ase()
    if correlation is None:
        raise RuntimeError("Correlation C++ extension module is not loaded.")

    is_periodic = any(atoms.pbc) or (hasattr(atoms, "cell") and getattr(atoms.cell, "volume", 0.0) > 1e-6)

    if is_periodic:
        cellpar = [float(p) for p in atoms.cell.cellpar()]
        cell = correlation.Cell(cellpar)
    else:
        cell = correlation.Cell()

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()

    for sym, pos in zip(symbols, positions):
        cell.add_atom(str(sym), [float(pos[0]), float(pos[1]), float(pos[2])])

    return cell


def to_ase(cell: Any) -> Any:
    """
    Convert a Correlation Cell into an ASE Atoms object.

    Parameters
    ----------
    cell : correlation.Cell
        Source Correlation simulation cell.

    Returns
    -------
    ase.Atoms
        Converted ASE Atoms object with matching positions, symbols, and cell parameters.

    Raises
    ------
    ImportError
        If ASE is not installed in the current environment.
    """
    Atoms = _require_ase()

    if hasattr(cell, "positions") and cell.positions is not None and len(cell.positions) > 0:
        positions = np.array(cell.positions, copy=True)
    else:
        positions = np.array([a.position for a in cell.atoms], dtype=float) if cell.atoms else np.empty((0, 3))

    symbols = [a.element.symbol for a in cell.atoms]
    lattice_params = cell.get_lattice_parameters() if hasattr(cell, "get_lattice_parameters") else None

    is_periodic = False
    cellpar = None
    if lattice_params is not None and len(lattice_params) == 6:
        if lattice_params[0] > 1e-6 and lattice_params[1] > 1e-6 and lattice_params[2] > 1e-6:
            cellpar = [float(p) for p in lattice_params]
            is_periodic = True

    return Atoms(symbols=symbols, positions=positions, cell=cellpar, pbc=is_periodic)


def from_ase_trajectory(images: Iterable[Any]) -> Any:
    """
    Convert an iterable or sequence of ASE Atoms frames into a Correlation Trajectory.

    Parameters
    ----------
    images : Iterable[ase.Atoms]
        Sequence of atomic snapshots (e.g., from ase.io.Trajectory).

    Returns
    -------
    correlation.Trajectory
        Constructed Trajectory containing all frames.
    """
    _require_ase()
    if correlation is None:
        raise RuntimeError("Correlation C++ extension module is not loaded.")

    traj = correlation.Trajectory()
    for img in images:
        frame = from_ase(img)
        if hasattr(traj, "add_frame"):
            traj.add_frame(frame)
        else:
            traj.frames.append(frame)
    return traj


def to_ase_trajectory(traj: Any) -> List[Any]:
    """
    Convert a Correlation Trajectory into a list of ASE Atoms objects.

    Parameters
    ----------
    traj : correlation.Trajectory
        Source multi-frame trajectory.

    Returns
    -------
    list[ase.Atoms]
        List of converted ASE Atoms snapshots.
    """
    _require_ase()
    num_frames = len(traj)
    return [to_ase(traj[i]) for i in range(num_frames)]


# -----------------------------------------------------------------------------
# Pymatgen Adapters
# -----------------------------------------------------------------------------

def from_pymatgen(structure_or_molecule: Any) -> Any:
    """
    Convert a Pymatgen Structure (periodic) or Molecule (non-periodic) into a Correlation Cell.

    Parameters
    ----------
    structure_or_molecule : pymatgen.core.Structure or pymatgen.core.Molecule
        Source Pymatgen representation.

    Returns
    -------
    correlation.Cell
        Constructed Correlation Cell.

    Raises
    ------
    ImportError
        If Pymatgen is not installed.
    """
    Structure, Molecule, _ = _require_pymatgen()
    if correlation is None:
        raise RuntimeError("Correlation C++ extension module is not loaded.")

    if isinstance(structure_or_molecule, Structure):
        lat = structure_or_molecule.lattice
        cell = correlation.Cell([float(lat.a), float(lat.b), float(lat.c),
                                 float(lat.alpha), float(lat.beta), float(lat.gamma)])
    else:
        cell = correlation.Cell()

    for site in structure_or_molecule:
        symbol = getattr(site.specie, "symbol", str(site.specie))
        coords = [float(site.coords[0]), float(site.coords[1]), float(site.coords[2])]
        cell.add_atom(symbol, coords)

    return cell


def to_pymatgen(cell: Any) -> Any:
    """
    Convert a Correlation Cell into a Pymatgen Structure (if periodic) or Molecule (if non-periodic).

    Parameters
    ----------
    cell : correlation.Cell
        Source Correlation simulation cell.

    Returns
    -------
    pymatgen.core.Structure or pymatgen.core.Molecule
        Converted Pymatgen data structure.

    Raises
    ------
    ImportError
        If Pymatgen is not installed.
    """
    Structure, Molecule, Lattice = _require_pymatgen()

    if hasattr(cell, "positions") and cell.positions is not None and len(cell.positions) > 0:
        coords = np.array(cell.positions, copy=True)
    else:
        coords = np.array([a.position for a in cell.atoms], dtype=float) if cell.atoms else np.empty((0, 3))

    symbols = [a.element.symbol for a in cell.atoms]
    params = cell.get_lattice_parameters() if hasattr(cell, "get_lattice_parameters") else None

    is_periodic = (
        params is not None
        and len(params) == 6
        and params[0] > 1e-6
        and params[1] > 1e-6
        and params[2] > 1e-6
    )

    if is_periodic:
        lattice = Lattice.from_parameters(
            float(params[0]), float(params[1]), float(params[2]),
            float(params[3]), float(params[4]), float(params[5])
        )
        return Structure(lattice, symbols, coords, coords_are_cartesian=True)

    return Molecule(symbols, coords)


# -----------------------------------------------------------------------------
# Dynamic Class Registration
# -----------------------------------------------------------------------------

def _register_adapters() -> None:
    """
    Attach adapter conversion methods directly onto correlation.Cell and correlation.Trajectory.
    """
    if correlation is None:
        return

    if hasattr(correlation, "Cell"):
        correlation.Cell.to_ase = to_ase
        correlation.Cell.from_ase = staticmethod(from_ase)
        correlation.Cell.to_pymatgen = to_pymatgen
        correlation.Cell.from_pymatgen = staticmethod(from_pymatgen)

    if hasattr(correlation, "Trajectory"):
        correlation.Trajectory.to_ase = to_ase_trajectory
        correlation.Trajectory.from_ase = staticmethod(from_ase_trajectory)
