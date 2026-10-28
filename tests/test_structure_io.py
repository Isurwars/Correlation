"""
test_structure_io.py — Tests for atomic structure readers, writers, and format converter bindings.
"""
import importlib
import os
import tempfile
import pytest

try:
    import _correlation as correlation
except ImportError:
    try:
        import matcorr as correlation
    except ImportError:
        correlation = None

if correlation is None or not hasattr(correlation, "read_structure"):
    pytest.skip("Correlation Python bindings not built", allow_module_level=True)


def test_read_structure_and_write(tmp_path):
    car_file = "tests/data/car/si_crystal.car"
    assert os.path.exists(car_file)

    cell = correlation.read_structure(car_file)
    assert len(cell) == 8
    assert cell.volume > 0.0

    out_vasp = tmp_path / "silicon.vasp"
    correlation.write(cell, str(out_vasp), format="poscar")
    assert out_vasp.exists()

    cell_from_vasp = correlation.read_structure(str(out_vasp))
    assert len(cell_from_vasp) == 8
    assert abs(cell_from_vasp.volume - cell.volume) < 1e-2


def test_cell_to_file_method(tmp_path):
    poscar_file = "tests/data/vasp/Si.poscar"
    assert os.path.exists(poscar_file)

    cell = correlation.read_structure(poscar_file)
    out_cell = tmp_path / "silicon.cell"
    cell.to_file(str(out_cell), format="cell")
    assert out_cell.exists()

    cell_from_cell = correlation.read_structure(str(out_cell))
    assert len(cell_from_cell) == 8


def test_convert_pipeline(tmp_path):
    car_file = "tests/data/car/si_crystal.car"
    out_cif = tmp_path / "silicon.cif"

    correlation.convert(car_file, str(out_cif))
    assert out_cif.exists()

    cell = correlation.read_structure(str(out_cif))
    assert len(cell) == 8
