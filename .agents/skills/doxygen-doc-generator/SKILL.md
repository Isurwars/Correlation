---
name: doxygen-doc-generator
description: Audits C++ headers for missing Doxygen documentation blocks, validates tag completeness (@param, @return, @throws), and formats compliant block comments.
---

# Doxygen Documentation Generator Skill

This skill enforces 100% Doxygen documentation compliance across all public/protected interfaces in C++ header files (`.hpp`, `.h`).

## 1. Compliance Standard

Every `class`, `struct`, `enum`, `concept`, and public/protected member function inside header files must declare a formatted Doxygen comment block.

### Format Specification
Use block style `/** ... */`. Annotate parameter directions explicitly:

```cpp
/**
 * @brief Computes the pair distribution function g(r) from atomistic positions.
 * 
 * Performs distance histogram binning over Cartesian particle coordinates
 * using thread-local OpenMP reduction buffers.
 * 
 * @param[in] coordinates Span of 3D Cartesian coordinates for all atoms in frame.
 * @param[in] cutoff Maximum radial distance limit in Angstroms (must be > 0).
 * @param[in] bin_width Radial bin resolution in Angstroms (must be > 0).
 * @return Histogram profile containing bin centers and g(r) amplitudes.
 * @throws std::invalid_argument If cutoff or bin_width is non-positive.
 */
[[nodiscard]] DistributionProfile compute_rdf(
    std::span<const Vector3D> coordinates,
    double cutoff,
    double bin_width
);
```

---

## 2. Tag Reference Matrix

| Doxygen Tag | Usage Requirement | Example |
| :--- | :--- | :--- |
| `@brief` | Required for all elements (1 sentence description) | `@brief Manages trajectory file parsing.` |
| `@param[in]` | Input parameter (read-only) | `@param[in] filename Path to input file.` |
| `@param[out]` | Output parameter (written by function) | `@param[out] result Target output buffer.` |
| `@param[in,out]` | Modified parameter | `@param[in,out] state Active calculator state.` |
| `@return` | Required for non-void returning functions | `@return Computed RDF histogram.` |
| `@throws` | Required if function can throw exceptions | `@throws std::runtime_error On parse failure.` |

---

## 3. Audit Protocol

When asked to audit or document a file:
1. Scan header file (`.hpp`) for public methods lacking `/**` blocks.
2. Verify all parameters match `@param` tags.
3. Verify return type matches `@return` tag.
4. Ensure no inline `NOLINT` or suppression tags are inserted.
