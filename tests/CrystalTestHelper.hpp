#pragma once

#include "core/Cell.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <string>

namespace correlation::testing::crystals {

/// Create a Face-Centered Cubic (FCC) supercell.
/// \param a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param nx Number of unit cells along X
/// \param ny Number of unit cells along Y
/// \param nz Number of unit cells along Z
inline correlation::core::Cell createFCCCell(real_t a, const std::string &element = "Cu", int nx = 2, int ny = 2,
                                             int nz = 2) {
  real_t const box_x = static_cast<real_t>(nx) * a;
  real_t const box_y = static_cast<real_t>(ny) * a;
  real_t const box_z = static_cast<real_t>(nz) * a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        real_t const ox = static_cast<real_t>(ix) * a;
        real_t const oy = static_cast<real_t>(iy) * a;
        real_t const oz = static_cast<real_t>(iz) * a;

        cell.addAtom(element, {ox, oy, oz});
        cell.addAtom(element, {ox + static_cast<real_t>(0.5) * a, oy + static_cast<real_t>(0.5) * a, oz});
        cell.addAtom(element, {ox + static_cast<real_t>(0.5) * a, oy, oz + static_cast<real_t>(0.5) * a});
        cell.addAtom(element, {ox, oy + static_cast<real_t>(0.5) * a, oz + static_cast<real_t>(0.5) * a});
      }
    }
  }

  return cell;
}

/// Create a Body-Centered Cubic (BCC) supercell.
/// \param a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param nx Number of unit cells along X
/// \param ny Number of unit cells along Y
/// \param nz Number of unit cells along Z
inline correlation::core::Cell createBCCCell(real_t a, const std::string &element = "Fe", int nx = 2, int ny = 2,
                                             int nz = 2) {
  real_t const box_x = static_cast<real_t>(nx) * a;
  real_t const box_y = static_cast<real_t>(ny) * a;
  real_t const box_z = static_cast<real_t>(nz) * a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        real_t const ox = static_cast<real_t>(ix) * a;
        real_t const oy = static_cast<real_t>(iy) * a;
        real_t const oz = static_cast<real_t>(iz) * a;

        cell.addAtom(element, {ox, oy, oz});
        cell.addAtom(element, {ox + static_cast<real_t>(0.5) * a, oy + static_cast<real_t>(0.5) * a,
                               oz + static_cast<real_t>(0.5) * a});
      }
    }
  }

  return cell;
}

/// Create a Diamond Cubic supercell.
/// \param a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param nx Number of unit cells along X
/// \param ny Number of unit cells along Y
/// \param nz Number of unit cells along Z
inline correlation::core::Cell createDiamondCell(real_t a, const std::string &element = "Si", int nx = 2, int ny = 2,
                                                 int nz = 2) {
  real_t const box_x = static_cast<real_t>(nx) * a;
  real_t const box_y = static_cast<real_t>(ny) * a;
  real_t const box_z = static_cast<real_t>(nz) * a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        real_t const ox = static_cast<real_t>(ix) * a;
        real_t const oy = static_cast<real_t>(iy) * a;
        real_t const oz = static_cast<real_t>(iz) * a;

        // FCC basis
        cell.addAtom(element, {ox, oy, oz});
        cell.addAtom(element, {ox + static_cast<real_t>(0.5) * a, oy + static_cast<real_t>(0.5) * a, oz});
        cell.addAtom(element, {ox + static_cast<real_t>(0.5) * a, oy, oz + static_cast<real_t>(0.5) * a});
        cell.addAtom(element, {ox, oy + static_cast<real_t>(0.5) * a, oz + static_cast<real_t>(0.5) * a});

        // Shifted FCC basis by (a/4, a/4, a/4)
        cell.addAtom(element, {ox + static_cast<real_t>(0.25) * a, oy + static_cast<real_t>(0.25) * a,
                               oz + static_cast<real_t>(0.25) * a});
        cell.addAtom(element, {ox + static_cast<real_t>(0.75) * a, oy + static_cast<real_t>(0.75) * a,
                               oz + static_cast<real_t>(0.25) * a});
        cell.addAtom(element, {ox + static_cast<real_t>(0.75) * a, oy + static_cast<real_t>(0.25) * a,
                               oz + static_cast<real_t>(0.75) * a});
        cell.addAtom(element, {ox + static_cast<real_t>(0.25) * a, oy + static_cast<real_t>(0.75) * a,
                               oz + static_cast<real_t>(0.75) * a});
      }
    }
  }

  return cell;
}

/// Create a Rock-Salt (NaCl) supercell.
/// \param a Lattice constant in Angstroms
/// \param cat_elem Cation element symbol (e.g. "Na")
/// \param an_elem Anion element symbol (e.g. "Cl")
/// \param nx Number of unit cells along X
/// \param ny Number of unit cells along Y
/// \param nz Number of unit cells along Z
inline correlation::core::Cell createNaClCell(real_t a, const std::string &cat_elem = "Na",
                                              const std::string &an_elem = "Cl", int nx = 2, int ny = 2, int nz = 2) {
  real_t const box_x = static_cast<real_t>(nx) * a;
  real_t const box_y = static_cast<real_t>(ny) * a;
  real_t const box_z = static_cast<real_t>(nz) * a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        real_t const ox = static_cast<real_t>(ix) * a;
        real_t const oy = static_cast<real_t>(iy) * a;
        real_t const oz = static_cast<real_t>(iz) * a;

        // Cations (Na) at FCC lattice sites
        cell.addAtom(cat_elem, {ox, oy, oz});
        cell.addAtom(cat_elem, {ox + static_cast<real_t>(0.5) * a, oy + static_cast<real_t>(0.5) * a, oz});
        cell.addAtom(cat_elem, {ox + static_cast<real_t>(0.5) * a, oy, oz + static_cast<real_t>(0.5) * a});
        cell.addAtom(cat_elem, {ox, oy + static_cast<real_t>(0.5) * a, oz + static_cast<real_t>(0.5) * a});

        // Anions (Cl) at FCC sites shifted by (a/2, 0, 0)
        cell.addAtom(an_elem, {ox + static_cast<real_t>(0.5) * a, oy, oz});
        cell.addAtom(an_elem, {ox, oy + static_cast<real_t>(0.5) * a, oz});
        cell.addAtom(an_elem, {ox, oy, oz + static_cast<real_t>(0.5) * a});
        cell.addAtom(an_elem, {ox + static_cast<real_t>(0.5) * a, oy + static_cast<real_t>(0.5) * a,
                               oz + static_cast<real_t>(0.5) * a});
      }
    }
  }

  return cell;
}

} // namespace correlation::testing::crystals
