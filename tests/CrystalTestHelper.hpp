#pragma once

#include "core/Cell.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"

#include <cmath>
#include <numbers>
#include <string>

namespace correlation::testing::crystals {

/// Create a Simple Cubic (SC) supercell.
/// \param a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param n_x Number of unit cells along X
/// \param n_y Number of unit cells along Y
/// \param n_z Number of unit cells along Z
inline correlation::core::Cell createSimpleCubicCell(real_t lat_a, const std::string &element = "Ar", int n_x = 2,
                                                     int n_y = 2, int n_z = 2) {
  real_t const box_x = static_cast<real_t>(n_x) * lat_a;
  real_t const box_y = static_cast<real_t>(n_y) * lat_a;
  real_t const box_z = static_cast<real_t>(n_z) * lat_a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int i_x = 0; i_x < n_x; ++i_x) {
    for (int i_y = 0; i_y < n_y; ++i_y) {
      for (int i_z = 0; i_z < n_z; ++i_z) {
        real_t const o_x = static_cast<real_t>(i_x) * lat_a;
        real_t const o_y = static_cast<real_t>(i_y) * lat_a;
        real_t const o_z = static_cast<real_t>(i_z) * lat_a;
        cell.addAtom(element, {o_x, o_y, o_z});
      }
    }
  }

  return cell;
}

/// Create a Face-Centered Cubic (FCC) supercell.
/// \param lat_a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param n_x Number of unit cells along X
/// \param n_y Number of unit cells along Y
/// \param n_z Number of unit cells along Z
inline correlation::core::Cell createFCCCell(real_t lat_a, const std::string &element = "Cu", int n_x = 2, int n_y = 2,
                                             int n_z = 2) {
  real_t const box_x = static_cast<real_t>(n_x) * lat_a;
  real_t const box_y = static_cast<real_t>(n_y) * lat_a;
  real_t const box_z = static_cast<real_t>(n_z) * lat_a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int i_x = 0; i_x < n_x; ++i_x) {
    for (int i_y = 0; i_y < n_y; ++i_y) {
      for (int i_z = 0; i_z < n_z; ++i_z) {
        real_t const o_x = static_cast<real_t>(i_x) * lat_a;
        real_t const o_y = static_cast<real_t>(i_y) * lat_a;
        real_t const o_z = static_cast<real_t>(i_z) * lat_a;

        cell.addAtom(element, {o_x, o_y, o_z});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.5) * lat_a, o_y + static_cast<real_t>(0.5) * lat_a, o_z});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.5) * lat_a, o_y, o_z + static_cast<real_t>(0.5) * lat_a});
        cell.addAtom(element, {o_x, o_y + static_cast<real_t>(0.5) * lat_a, o_z + static_cast<real_t>(0.5) * lat_a});
      }
    }
  }

  return cell;
}

/// Create a Body-Centered Cubic (BCC) supercell.
/// \param lat_a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param n_x Number of unit cells along X
/// \param n_y Number of unit cells along Y
/// \param n_z Number of unit cells along Z
inline correlation::core::Cell createBCCCell(real_t lat_a, const std::string &element = "Fe", int n_x = 2, int n_y = 2,
                                             int n_z = 2) {
  real_t const box_x = static_cast<real_t>(n_x) * lat_a;
  real_t const box_y = static_cast<real_t>(n_y) * lat_a;
  real_t const box_z = static_cast<real_t>(n_z) * lat_a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int i_x = 0; i_x < n_x; ++i_x) {
    for (int i_y = 0; i_y < n_y; ++i_y) {
      for (int i_z = 0; i_z < n_z; ++i_z) {
        real_t const o_x = static_cast<real_t>(i_x) * lat_a;
        real_t const o_y = static_cast<real_t>(i_y) * lat_a;
        real_t const o_z = static_cast<real_t>(i_z) * lat_a;

        cell.addAtom(element, {o_x, o_y, o_z});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.5) * lat_a, o_y + static_cast<real_t>(0.5) * lat_a,
                               o_z + static_cast<real_t>(0.5) * lat_a});
      }
    }
  }

  return cell;
}

/// Create a Hexagonal Close-Packed (HCP) supercell.
/// \param lat_a Lattice constant in basal plane in Angstroms
/// \param lat_c Lattice constant along Z axis in Angstroms (typically c/a ~ 1.633)
/// \param element Chemical element symbol
/// \param n_x Number of unit cells along X
/// \param n_y Number of unit cells along Y
/// \param n_z Number of unit cells along Z
inline correlation::core::Cell createHCPCell(real_t lat_a, real_t lat_c, const std::string &element = "Mg", int n_x = 2,
                                             int n_y = 2, int n_z = 2) {
  constexpr real_t SQRT3 = std::numbers::sqrt3_v<real_t>;
  real_t const box_x = static_cast<real_t>(n_x) * lat_a;
  real_t const box_y = static_cast<real_t>(n_y) * lat_a * SQRT3 / static_cast<real_t>(2.0);
  real_t const box_z = static_cast<real_t>(n_z) * lat_c;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int i_x = 0; i_x < n_x; ++i_x) {
    for (int i_y = 0; i_y < n_y; ++i_y) {
      for (int i_z = 0; i_z < n_z; ++i_z) {
        real_t const o_x = static_cast<real_t>(i_x) * lat_a +
                           (i_y % 2 == 1 ? static_cast<real_t>(0.5) * lat_a : static_cast<real_t>(0.0));
        real_t const o_y = static_cast<real_t>(i_y) * lat_a * SQRT3 / static_cast<real_t>(2.0);
        real_t const o_z = static_cast<real_t>(i_z) * lat_c;

        cell.addAtom(element, {o_x, o_y, o_z});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.5) * lat_a, o_y + lat_a / static_cast<real_t>(2.0 * SQRT3),
                               o_z + static_cast<real_t>(0.5) * lat_c});
      }
    }
  }

  return cell;
}

/// Create a Diamond Cubic supercell.
/// \param lat_a Lattice constant in Angstroms
/// \param element Chemical element symbol
/// \param n_x Number of unit cells along X
/// \param n_y Number of unit cells along Y
/// \param n_z Number of unit cells along Z
inline correlation::core::Cell createDiamondCell(real_t lat_a, const std::string &element = "Si", int n_x = 2,
                                                 int n_y = 2, int n_z = 2) {
  real_t const box_x = static_cast<real_t>(n_x) * lat_a;
  real_t const box_y = static_cast<real_t>(n_y) * lat_a;
  real_t const box_z = static_cast<real_t>(n_z) * lat_a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int i_x = 0; i_x < n_x; ++i_x) {
    for (int i_y = 0; i_y < n_y; ++i_y) {
      for (int i_z = 0; i_z < n_z; ++i_z) {
        real_t const o_x = static_cast<real_t>(i_x) * lat_a;
        real_t const o_y = static_cast<real_t>(i_y) * lat_a;
        real_t const o_z = static_cast<real_t>(i_z) * lat_a;

        // FCC basis
        cell.addAtom(element, {o_x, o_y, o_z});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.5) * lat_a, o_y + static_cast<real_t>(0.5) * lat_a, o_z});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.5) * lat_a, o_y, o_z + static_cast<real_t>(0.5) * lat_a});
        cell.addAtom(element, {o_x, o_y + static_cast<real_t>(0.5) * lat_a, o_z + static_cast<real_t>(0.5) * lat_a});

        // Shifted FCC basis by (a/4, a/4, a/4)
        cell.addAtom(element, {o_x + static_cast<real_t>(0.25) * lat_a, o_y + static_cast<real_t>(0.25) * lat_a,
                               o_z + static_cast<real_t>(0.25) * lat_a});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.75) * lat_a, o_y + static_cast<real_t>(0.75) * lat_a,
                               o_z + static_cast<real_t>(0.25) * lat_a});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.75) * lat_a, o_y + static_cast<real_t>(0.25) * lat_a,
                               o_z + static_cast<real_t>(0.75) * lat_a});
        cell.addAtom(element, {o_x + static_cast<real_t>(0.25) * lat_a, o_y + static_cast<real_t>(0.75) * lat_a,
                               o_z + static_cast<real_t>(0.75) * lat_a});
      }
    }
  }

  return cell;
}

/// Create a Rock-Salt (NaCl) supercell.
/// \param a Lattice constant in Angstroms
/// \param cat_elem Cation element symbol (e.g. "Na")
/// \param an_elem Anion element symbol (e.g. "Cl")
/// \param n_x Number of unit cells along X
/// \param n_y Number of unit cells along Y
/// \param n_z Number of unit cells along Z
inline correlation::core::Cell createNaClCell(real_t lat_a, const std::string &cat_elem = "Na",
                                              const std::string &an_elem = "Cl", int n_x = 2, int n_y = 2,
                                              int n_z = 2) {
  real_t const box_x = static_cast<real_t>(n_x) * lat_a;
  real_t const box_y = static_cast<real_t>(n_y) * lat_a;
  real_t const box_z = static_cast<real_t>(n_z) * lat_a;

  correlation::core::Cell cell({box_x, box_y, box_z, 90.0, 90.0, 90.0});

  for (int i_x = 0; i_x < n_x; ++i_x) {
    for (int i_y = 0; i_y < n_y; ++i_y) {
      for (int i_z = 0; i_z < n_z; ++i_z) {
        real_t const o_x = static_cast<real_t>(i_x) * lat_a;
        real_t const o_y = static_cast<real_t>(i_y) * lat_a;
        real_t const o_z = static_cast<real_t>(i_z) * lat_a;

        // Cations (Na) at FCC lattice sites
        cell.addAtom(cat_elem, {o_x, o_y, o_z});
        cell.addAtom(cat_elem, {o_x + static_cast<real_t>(0.5) * lat_a, o_y + static_cast<real_t>(0.5) * lat_a, o_z});
        cell.addAtom(cat_elem, {o_x + static_cast<real_t>(0.5) * lat_a, o_y, o_z + static_cast<real_t>(0.5) * lat_a});
        cell.addAtom(cat_elem, {o_x, o_y + static_cast<real_t>(0.5) * lat_a, o_z + static_cast<real_t>(0.5) * lat_a});

        // Anions (Cl) at FCC sites shifted by (a/2, 0, 0)
        cell.addAtom(an_elem, {o_x + static_cast<real_t>(0.5) * lat_a, o_y, o_z});
        cell.addAtom(an_elem, {o_x, o_y + static_cast<real_t>(0.5) * lat_a, o_z});
        cell.addAtom(an_elem, {o_x, o_y, o_z + static_cast<real_t>(0.5) * lat_a});
        cell.addAtom(an_elem, {o_x + static_cast<real_t>(0.5) * lat_a, o_y + static_cast<real_t>(0.5) * lat_a,
                               o_z + static_cast<real_t>(0.5) * lat_a});
      }
    }
  }

  return cell;
}

/// Parameters for 2-atom Dimer cell creation.
struct DimerCellOptions {
  std::string elem1{"Ar"};
  std::string elem2{"Ar"};
  real_t dist{1.5};
  real_t box_size{10.0};
};

/// Create a 2-atom Dimer cell.
/// \param opts Configuration options for the dimer cell
inline correlation::core::Cell createDimerCell(DimerCellOptions const &opts = {}) {
  correlation::core::Cell cell({opts.box_size, opts.box_size, opts.box_size, 90.0, 90.0, 90.0});
  real_t const mid = opts.box_size * static_cast<real_t>(0.5);
  cell.addAtom(opts.elem1, {mid, mid, mid});
  cell.addAtom(opts.elem2, {mid + opts.dist, mid, mid});
  return cell;
}

/// Parameters for Water (H2O) molecule cell creation.
struct WaterMoleculeOptions {
  correlation::math::Vector3<real_t> O_pos{5.0, 5.0, 5.0};
  real_t r_OH{0.957};
  real_t angle_HOH_deg{104.52};
  real_t box_size{10.0};
};

/// Create a Water (H2O) molecule cell.
/// \param opts Configuration options for the water molecule cell
inline correlation::core::Cell createWaterMoleculeCell(WaterMoleculeOptions const &opts = {}) {
  correlation::core::Cell cell({opts.box_size, opts.box_size, opts.box_size, 90.0, 90.0, 90.0});
  cell.addAtom("O", opts.O_pos);

  real_t const ang_rad = opts.angle_HOH_deg * static_cast<real_t>(correlation::math::deg_to_rad);
  // H1 along +X axis from O
  cell.addAtom("H", {opts.O_pos.x() + opts.r_OH, opts.O_pos.y(), opts.O_pos.z()});
  // H2 rotated by ang_rad in XY plane
  cell.addAtom("H", {opts.O_pos.x() + opts.r_OH * std::cos(ang_rad), opts.O_pos.y() + opts.r_OH * std::sin(ang_rad),
                     opts.O_pos.z()});

  return cell;
}

/// Parameters for Triatomic Bent / Angle molecule cell creation.
struct TriatomicAngleCellOptions {
  std::string center_elem{"O"};
  std::string arm_elem1{"H"};
  std::string arm_elem2{"H"};
  real_t dist1{0.957};
  real_t dist2{0.957};
  real_t angle_deg{104.52};
  real_t box_size{20.0};
};

/// Create a Triatomic Bent / Angle molecule cell.
/// \param opts Configuration options for the triatomic angle cell
inline correlation::core::Cell createTriatomicAngleCell(TriatomicAngleCellOptions const &opts = {}) {
  correlation::core::Cell cell({opts.box_size, opts.box_size, opts.box_size, 90.0, 90.0, 90.0});
  real_t const mid = opts.box_size * static_cast<real_t>(0.5);
  cell.addAtom(opts.center_elem, {mid, mid, mid});
  cell.addAtom(opts.arm_elem1, {mid + opts.dist1, mid, mid});

  real_t const ang_rad = opts.angle_deg * static_cast<real_t>(correlation::math::deg_to_rad);
  cell.addAtom(opts.arm_elem2, {mid + opts.dist2 * std::cos(ang_rad), mid + opts.dist2 * std::sin(ang_rad), mid});
  return cell;
}


/// Parameters for Regular N-Polygon Ring cell creation.
struct RingCellOptions {
  size_t num_atoms{6};
  std::string element{"C"};
  real_t radius{1.4};
  real_t box_size{20.0};
};

/// Create a Regular N-Polygon Ring cell.
/// \param opts Configuration options for the ring cell
inline correlation::core::Cell createRingCell(RingCellOptions const &opts = {}) {
  correlation::core::Cell cell({opts.box_size, opts.box_size, opts.box_size, 90.0, 90.0, 90.0});
  real_t const mid = opts.box_size * static_cast<real_t>(0.5);

  real_t const d_theta = static_cast<real_t>(correlation::math::two_pi) / static_cast<real_t>(opts.num_atoms);
  for (size_t i = 0; i < opts.num_atoms; ++i) {
    real_t const theta = static_cast<real_t>(i) * d_theta;
    cell.addAtom(opts.element, {mid + opts.radius * std::cos(theta), mid + opts.radius * std::sin(theta), mid});
  }

  return cell;
}

/// Parameters for 13-atom Icosahedral Cluster creation.
struct IcosahedralClusterOptions {
  std::string center_elem{"Ar"};
  std::string shell_elem{"Ar"};
  real_t r_bond{1.0};
  real_t box_size{10.0};
};

/// Create a 13-atom Icosahedral Cluster cell.
/// \param opts Configuration options for the cluster cell
inline correlation::core::Cell createIcosahedralClusterCell(IcosahedralClusterOptions const &opts = {}) {
  correlation::core::Cell cell({opts.box_size, opts.box_size, opts.box_size, 90.0, 90.0, 90.0});
  real_t const mid = opts.box_size * static_cast<real_t>(0.5);

  // Central atom
  cell.addAtom(opts.center_elem, {mid, mid, mid});

  real_t const phi = std::numbers::phi_v<real_t>;
  real_t const scale = opts.r_bond / std::sqrt(static_cast<real_t>(1.0) + phi * phi);
  real_t const phi_scale = phi * scale;

  // 12 shell vertices
  cell.addAtom(opts.shell_elem, {mid, mid + scale, mid + phi_scale});
  cell.addAtom(opts.shell_elem, {mid, mid + scale, mid - phi_scale});
  cell.addAtom(opts.shell_elem, {mid, mid - scale, mid + phi_scale});
  cell.addAtom(opts.shell_elem, {mid, mid - scale, mid - phi_scale});

  cell.addAtom(opts.shell_elem, {mid + scale, mid + phi_scale, mid});
  cell.addAtom(opts.shell_elem, {mid + scale, mid - phi_scale, mid});
  cell.addAtom(opts.shell_elem, {mid - scale, mid + phi_scale, mid});
  cell.addAtom(opts.shell_elem, {mid - scale, mid - phi_scale, mid});

  cell.addAtom(opts.shell_elem, {mid + phi_scale, mid, mid + scale});
  cell.addAtom(opts.shell_elem, {mid + phi_scale, mid, mid - scale});
  cell.addAtom(opts.shell_elem, {mid - phi_scale, mid, mid + scale});
  cell.addAtom(opts.shell_elem, {mid - phi_scale, mid, mid - scale});

  return cell;
}

} // namespace correlation::testing::crystals
