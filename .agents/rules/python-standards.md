# Rule: Python Standards, Type Checking & Bindings Interop

*Activation Mode: Glob (`**/*.py`, `**/*.pyi`, `pyproject.toml`)*

## 1. Type Checking & Language Server Configuration
- **Type Checker:** Pyright / BasedPyright with `typeCheckingMode = "standard"`.
- **C-Extension Stubs:** The C++ module (`_correlation`) is a compiled shared library. All public interfaces must declare accurate typing stubs in `python/correlation/_correlation.pyi`.
- **Module Resolution:** Keep `reportMissingModuleSource = false` in `[tool.pyright]` and `[tool.basedpyright]` within `pyproject.toml` to suppress warnings on compiled extensions with stubs.
- **Inline Ignore Directive:** When re-exporting the compiled extension in `python/correlation/__init__.py`, use `# pyright: ignore[reportMissingModuleSource]`.

## 2. Python Testing Standards
- **Framework:** `pytest` is the standard Python test runner.
- **Test Discovery:** Execute tests inside the virtual environment:
  ```bash
  .venv/bin/pytest tests/test_*.py
  ```
- **NumPy Zero-Copy Invariants:** Test that array buffers passed between Python and C++ avoid redundant copying and verify dimensions/strides.

## 3. Package Management & Installation
- **Build Backend:** `scikit-build-core` with `pybind11`.
- **Development Installation:** Use editable / local build mode when testing bindings:
  ```bash
  pip install --no-build-isolation -e .
  ```
- **Python Version Compatibility:** Python 3.10 through 3.14.
