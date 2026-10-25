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
        __version__ = "4.1.0"

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
    "compute_ten_wolde_frenkel",
    "compute_tetrahedral_order",
    "compute_local_entropy",
    "compute_g6",
    "fit_correlation_length",
    "TenWoldeFrenkelParams",
    "TetrahedralParams",
    "LocalEntropyParams",
    "CorrelationLengthFit",
    "adapters",
    "__version__",
]
