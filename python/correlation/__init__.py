"""
Correlation — High-performance structural analysis for atomistic simulations.

This package provides Python bindings for the Correlation C++ analysis engine,
enabling computation of radial distribution functions (RDF), structure factors
S(Q), bond angle distributions, and other structural properties of liquids,
amorphous solids, and crystalline materials.

Example usage::

    import correlation

    # Read a structure file
    cell = correlation.Cell.from_file("structure.car")

    # Compute distribution functions
    df = correlation.DistributionFunctions(cell, cutoff=10.0)
    df.calculate_rdf(r_max=20.0, bin_width=0.05)

    # Access results
    hist = df.get_histogram("g_r")
    print(hist.bins, hist.partials)
"""

try:
    from correlation._correlation import *  # noqa: F401,F403
except ImportError as e:
    raise ImportError(
        "Failed to import the Correlation C++ extension module. "
        "Make sure the package was built correctly with: pip install ."
    ) from e
