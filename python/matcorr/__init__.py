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
        __version__ = "4.0.2"

# Import all core capabilities from correlation engine
try:
    import correlation
    from correlation import *  # noqa: F401, F403
    try:
        from correlation import adapters  # noqa: F401
    except ImportError:
        adapters = None
except ImportError:
    adapters = None

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
    "adapters",
    "__version__",
]
