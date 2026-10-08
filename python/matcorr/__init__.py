"""
matcorr — Materials Correlation Analysis Suite.

High-performance structural and correlation analysis for liquid, amorphous,
and nanostructured materials.
"""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("matcorr")
except PackageNotFoundError:
    try:
        __version__ = version("correlation-analysis")
    except PackageNotFoundError:
        __version__ = "3.9.11"

# Import all core capabilities from correlation engine
import correlation
from correlation import *  # noqa: F401, F403

__all__ = [
    "Cell",
    "Atom",
    "Trajectory",
    "TrajectoryAnalyzer",
    "StructureAnalyzer",
    "DistributionFunctions",
    "AnalysisSettings",
    "read",
    "write_csv",
    "get_calculator",
    "list_calculators",
    "get_writer",
    "list_writers",
    "__version__",
]
