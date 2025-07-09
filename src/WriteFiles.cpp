// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
#include "../include/WriteFiles.hpp"

#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

//---------------------------------------------------------------------------//
//------------------------------ Core Implementation ------------------------//
//---------------------------------------------------------------------------//

// Helper structure for output configuration
struct OutputConfig {
  std::string_view name;
  std::string_view suffix;
  std::function<std::vector<std::vector<double>>(const DistributionFunctions &)>
      get_data;
  std::function<std::string(const Cell &)> header_gen;
};

// Header generators
std::string pair_header(const Cell &cell, std::string unit,
                        std::string quantity) {
  std::stringstream ss;
  const int n = cell.elements().size();
  ss << std::right << std::setw(14) << "r (Å),";
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      ss << std::right << std::setw(14)
         << cell.elements()[i] + "-" + cell.elements()[j] + " " + unit << ",";
    }
  }
  ss << std::right << std::setw(12) << quantity;
  return ss.str();
}

std::string angle_header(const Cell &cell, std::string unit,
                         std::string quantity) {
  std::stringstream ss;
  const int n = cell.elements().size();
  ss << std::right << std::setw(14) << "theta (deg),";
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = j; k < n; ++k) {
        ss << std::right << std::setw(13)
           << cell.elements()[j] + "-" + cell.elements()[i] + "-" +
                  cell.elements()[k] + " " + unit
           << ",";
      }
    }
  }
  ss << std::right << std::setw(12) << quantity;
  return ss.str();
}

std::string coordination_header(const Cell &cell) {
  std::stringstream ss;
  const int n = cell.elements().size();
  ss << std::right << std::setw(14) << "Counts (#),";
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ss << std::right << std::setw(14)
         << cell.elements()[i] + " by " + cell.elements()[j] + ",";
    }
  }
  for (int i = 0; i < n; i++) {
    ss << std::right << std::setw(13) << cell.elements()[i] << " by any";
    if (i < n - 1)
      ss << ",";
  }

  return ss.str();
}

// Output configurations
const std::array<OutputConfig, 7> base_configs = {{
    {"J(r)", "_J.csv", [](auto &df) { return df.J(); },
     [](auto &cell) { return pair_header(cell, "(Å⁻¹)", "J(r) (Å⁻¹)"); }},

    {"g(r)", "_g.csv", [](auto &df) { return df.g(); },
     [](auto &cell) { return pair_header(cell, "(Å⁻¹)", "g(r) (Å⁻¹)"); }},

    {"G(r)", "_G.csv", [](auto &df) { return df.G(); },
     [](auto &cell) { return pair_header(cell, "(Å⁻¹)", "G(r) (Å⁻¹)"); }},

    {"F(theta)", "_PAD.csv", [](auto &df) { return df.F(); },
     [](auto &cell) {
       return angle_header(cell, "(deg⁻¹)", "F(theta) (deg⁻¹)");
     }},

    {"S(Q)", "_S.csv", [](auto &df) { return df.S(); },
     [](auto &cell) { return pair_header(cell, "(Å)", "S(q) (Å)"); }},

    {"XRD", "_XRD.csv", [](auto &df) { return df.X(); },
     [](auto &cell) { return pair_header(cell, "(deg)", "XRD (deg⁻¹)"); }},

    {"Coordination", "_Z.csv", [](auto &df) { return df.Z(); },
     coordination_header},
}};

const std::array<OutputConfig, 5> smoothed_configs = {
    {{"J(r) smoothed", "_J_smoothed.csv",
      [](auto &df) { return df.J_smoothed(); },
      [](auto &cell) { return pair_header(cell, "(Å⁻¹)", "J(r) (Å⁻¹)"); }},

     {"g(r) smoothed", "_g_smoothed.csv",
      [](auto &df) { return df.g_smoothed(); },
      [](auto &cell) { return pair_header(cell, "(Å⁻¹)", "g(r) (Å⁻¹)"); }},

     {"G(r) smoothed", "_G_smoothed.csv",
      [](auto &df) { return df.G_smoothed(); },
      [](auto &cell) { return pair_header(cell, "(Å⁻¹)", "G(r) (Å⁻¹)"); }},

     {"F(theta) smoothed", "_PAD_smoothed.csv",
      [](auto &df) { return df.F_smoothed(); },
      [](auto &cell) {
        return angle_header(cell, "(deg⁻¹)", "F(theta) (deg⁻¹)");
      }},

     {"S(Q) smoothed", "_S_smoothed.csv",
      [](auto &df) { return df.S_smoothed(); },
      [](auto &cell) { return pair_header(cell, "(Å)", "S(q) Å"); }}}};

template <typename T>
void write_impl(const std::vector<std::vector<T>> &data,
                const std::string &filename, const std::string &header) {
  std::ofstream file(filename, std::ios::binary);
  if (!file) {
    throw std::runtime_error("Failed to open file: " + filename);
  }

  // Write UTF-8 BOM
  const unsigned char bom[] = {0xEF, 0xBB, 0xBF};
  file.write(reinterpret_cast<const char *>(bom), sizeof(bom));
  file << std::fixed << std::setprecision(5);
  file << header << '\n';

  for (size_t j = 0; j < data[0].size(); ++j) {
    for (size_t i = 0; i < data.size(); ++i) {
      file << std::setw(12) << std::right << data[i][j];
      if (i < data.size() - 1)
        file << ",";
    }
    file << '\n';
  }

  if (!file) {
    throw std::runtime_error("Error writing to file: " + filename);
  }
}

//---------------------------------------------------------------------------//
//------------------------------- Write CSV ---------------------------------//
//---------------------------------------------------------------------------//

void WriteCSV(const DistributionFunctions &df, const std::string &base_path,
              bool smoothed) {
  const auto &cell = df.cell();
  const int n = cell.elements().size();

  // Write base distributions
  for (const auto &config : base_configs) {
    try {
      const auto data = config.get_data(df);
      const std::string filename = base_path + std::string(config.suffix);
      write_impl(data, filename, config.header_gen(cell));
      std::cout << "Writing " << config.name << '\n';
    } catch (const std::exception &e) {
      std::cerr << "Error writing " << config.name << ": " << e.what() << '\n';
    }
  }

  // Write smoothed versions if requested
  if (smoothed) {
    for (const auto &config : smoothed_configs) {
      try {
        const auto data = config.get_data(df);
        const std::string filename = base_path + std::string(config.suffix);
        write_impl(data, filename, config.header_gen(cell));
        std::cout << "Writing " << config.name << '\n';
      } catch (const std::exception &e) {
        std::cerr << "Error writing " << config.name << ": " << e.what()
                  << '\n';
      }
    }
  }
}
