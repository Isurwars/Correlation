// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <map>
#include <string>

namespace Correlation {
namespace WriterUtils {

// Define metadata structure
struct FunctionMetadata {
  std::string bin_label;
  std::string bin_unit;
  std::string data_unit;
  std::string description;
};

// Metadata mapping
inline const std::map<std::string, FunctionMetadata> metadata_map = {
    {"g_r", {"r (Å)", "Å", "Å^-1", "Radial Distribution Function"}},
    {"J_r", {"r (Å)", "Å", "Å^-1", "Radial Distribution of Electron Density"}},
    {"G_r", {"r (Å)", "Å", "Å^-1", "Reduced Radial Distribution Function"}},
    {"BAD",
     {"Bond Angle (°)", "Degrees", "degree^-1", "Bond Angle Distribution"}},
    {"DAD",
     {"Dihedral Angle (°)", "Degrees", "degree^-1",
      "Dihedral Angle Distribution"}},
    {"RD", {"Ring Size", "atoms", "counts", "Ring Distribution"}},
    {"S_q", {"q (Å^-1)", "Å^-1", "arbitrary units", "Structure Factor S(Q)"}},
    {"XRD",
     {"2theta", "Degrees (2theta)", "Intensity", "X-Ray Diffraction Pattern"}},
    {"CN", {"counts", "neighbors", "Count", "Coordination Number"}},
    {"VACF",
     {"Time (fs)", "fs", "Å^2/fs^2", "Velocity Autocorrelation Function"}},
    {"Normalized VACF",
     {"Time (fs)", "fs", "normalized",
      "Normalized Velocity Autocorrelation Function"}},
    {"VDOS",
     {"Frequency_THz", "THz", "arbitrary units",
      "Vibrational Density of States"}},
    {"Frequency_THz", {"Frequency_THz", "THz", "THz", "Frequency_THz"}},
    {"Frequency_cm_1",
     {"Frequency_cm_1", "cm^-1", "arbitrary units",
      "Frequency in wavenumbers"}},
    {"Frequency_meV",
     {"Frequency_meV", "meV", "arbitrary units", "Frequency in meV"}},
    {"MSD", {"Time (fs)", "fs", "Å²", "Mean Squared Displacement"}},
    {"D_eff",
     {"Time (fs)", "fs", "Å²/fs",
      "Running Diffusion Coefficient D(t) = MSD(t) / (6t)"}}};

} // namespace WriterUtils
} // namespace Correlation
