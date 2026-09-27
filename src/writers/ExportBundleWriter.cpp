/**
 * @file ExportBundleWriter.cpp
 * @brief Implementation of the consolidated ZIP export bundle writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/ExportBundleWriter.hpp"
#include "plotters/SvgHistogramRenderer.hpp"
#include "utils/ZipArchiveWriter.hpp"
#include "writers/CSVWriter.hpp"

#ifdef CORRELATION_USE_HDF5
#include "writers/HDF5Writer.hpp"
#endif

#ifdef CORRELATION_USE_ARROW
#include "writers/ArrowWriter.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <format>
#include <fstream>
#include <sstream>

namespace correlation::writers {

namespace {

bool isAlgorithmSelected(const std::string &hist_name, const std::vector<std::string> &selected) {
  if (selected.empty()) {
    return true;
  }
  std::string base_name = hist_name;
  if (base_name.ends_with("_raw")) {
    base_name = base_name.substr(0, base_name.length() - 4);
  }
  return std::ranges::find(selected, base_name) != selected.end() ||
         std::ranges::find(selected, hist_name) != selected.end();
}

std::string generateSummaryText(const correlation::analysis::DistributionFunctions &dists) {
  std::ostringstream summary_stream;
  summary_stream << "Correlation Analysis Summary\n";
  summary_stream << "============================\n\n";
  summary_stream << "Calculated Dynamic Properties:\n";
  summary_stream << "------------------------------\n";

  const real_t d_msd = dists.getDiffusionCoefficientMSD();
  const real_t d_vacf = dists.getDiffusionCoefficientVACF();
  const real_t tau = dists.getRelaxationTime();
  const real_t deb = dists.getDeborahNumber();

  summary_stream << "Self-diffusion coefficient (from MSD): ";
  if (d_msd > 0.0) {
    summary_stream << d_msd << " Å²/fs\n";
  } else {
    summary_stream << "N/A\n";
  }

  summary_stream << "Self-diffusion coefficient (from VACF): ";
  if (d_vacf > 0.0) {
    summary_stream << d_vacf << " Å²/fs\n";
  } else {
    summary_stream << "N/A\n";
  }

  summary_stream << "Relaxation time (from VACF): ";
  if (tau > 0.0) {
    summary_stream << tau << " fs\n";
  } else {
    summary_stream << "N/A\n";
  }

  summary_stream << "Deborah number: ";
  if (deb > 0.0) {
    summary_stream << deb << "\n";
  } else {
    summary_stream << "N/A\n";
  }

  return summary_stream.str();
}

std::expected<void, std::string>
packagePlots(utils::ZipArchiveWriter &zip,
             const correlation::analysis::DistributionFunctions &dists,
             const ExportBundleOptions &options) {
  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
      continue;
    }
    if (!isAlgorithmSelected(name, options.selected_algorithms)) {
      continue;
    }
    // Render in-memory SVG
    const std::string svg_str =
        plotters::renderHistogramAsSvg(hist, options.plot_config, {}, dists.getAshcroftWeights());

    std::string filename = hist.file_suffix;
    if (filename.starts_with('_')) {
      filename = filename.substr(1);
    }
    const std::string entry_path = std::format("plots/{}.svg", filename);
    auto res = zip.addFileFromString(entry_path, svg_str);
    if (!res) {
      return res;
    }
  }
  return {};
}

std::expected<void, std::string>
packageCSVs(utils::ZipArchiveWriter &zip, const correlation::analysis::DistributionFunctions &dists,
            const ExportBundleOptions &options) {
  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
      continue;
    }
    if (!isAlgorithmSelected(name, options.selected_algorithms)) {
      continue;
    }

    const correlation::analysis::Histogram *raw_companion = nullptr;
    const std::string raw_key = name + "_raw";
    if ((name == "PAD" || name == "DAD") && all_histograms.contains(raw_key)) {
      raw_companion = &all_histograms.at(raw_key);
    }

    std::ostringstream csv_ss;
    CSVWriter::writeHistogramToStream(csv_ss, hist, raw_companion);

    std::string filename = hist.file_suffix;
    if (filename.starts_with('_')) {
      filename = filename.substr(1);
    }
    const std::string entry_path = std::format("csv/{}.csv", filename);
    auto res = zip.addFileFromString(entry_path, csv_ss.str());
    if (!res) {
      return res;
    }
  }
  return {};
}

#ifdef CORRELATION_USE_HDF5
std::expected<void, std::string>
packageHDF5(utils::ZipArchiveWriter &zip, const correlation::analysis::DistributionFunctions &dists,
            const std::filesystem::path &temp_dir) {
  const auto temp_h5 = temp_dir / "analysis.h5";
  HDF5Writer::writeHDF(temp_h5.string(), dists);
  auto res = zip.addFileFromDisk("hdf5/analysis.h5", temp_h5);
  std::error_code ec;
  std::filesystem::remove(temp_h5, ec);
  return res;
}
#endif

#ifdef CORRELATION_USE_ARROW
std::expected<void, std::string>
packageParquet(utils::ZipArchiveWriter &zip,
               const correlation::analysis::DistributionFunctions &dists,
               const ExportBundleOptions &options, const std::filesystem::path &temp_dir) {
  const auto base_temp = temp_dir / "parquet_export";
  ArrowWriter::writeAllParquet(base_temp.string(), dists, options.smoothing);

  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    if (hist.file_suffix.empty()) {
      continue;
    }
    if (!isAlgorithmSelected(name, options.selected_algorithms)) {
      continue;
    }
    const auto src_parquet = temp_dir / std::format("parquet_export{}.parquet", hist.file_suffix);
    if (std::filesystem::exists(src_parquet)) {
      std::string filename = hist.file_suffix;
      if (filename.starts_with('_')) {
        filename = filename.substr(1);
      }
      const std::string entry_path = std::format("parquet/{}.parquet", filename);
      auto res = zip.addFileFromDisk(entry_path, src_parquet);
      std::error_code ec;
      std::filesystem::remove(src_parquet, ec);
      if (!res) {
        return res;
      }
    }
  }
  return {};
}
#endif

void writeFolderPlots(const std::filesystem::path &plots_dir,
                      const correlation::analysis::DistributionFunctions &dists,
                      const ExportBundleOptions &options) {
  std::error_code dir_error;
  std::filesystem::create_directories(plots_dir, dir_error);
  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
      continue;
    }
    if (!isAlgorithmSelected(name, options.selected_algorithms)) {
      continue;
    }
    const std::string svg_str =
        plotters::renderHistogramAsSvg(hist, options.plot_config, {}, dists.getAshcroftWeights());
    std::string filename = hist.file_suffix;
    if (filename.starts_with('_')) {
      filename = filename.substr(1);
    }
    const auto out_path = plots_dir / std::format("{}.svg", filename);
    std::ofstream out(out_path);
    if (out) {
      out << svg_str;
    }
  }
}

void writeFolderCSVs(const std::filesystem::path &csv_dir,
                     const correlation::analysis::DistributionFunctions &dists,
                     const ExportBundleOptions &options) {
  std::error_code dir_error;
  std::filesystem::create_directories(csv_dir, dir_error);
  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    if (hist.partials.empty() || hist.bins.empty() || hist.file_suffix.empty()) {
      continue;
    }
    if (!isAlgorithmSelected(name, options.selected_algorithms)) {
      continue;
    }
    const correlation::analysis::Histogram *raw_companion = nullptr;
    const std::string raw_key = name + "_raw";
    if ((name == "PAD" || name == "DAD") && all_histograms.contains(raw_key)) {
      raw_companion = &all_histograms.at(raw_key);
    }
    std::string filename = hist.file_suffix;
    if (filename.starts_with('_')) {
      filename = filename.substr(1);
    }
    const auto out_path = csv_dir / std::format("{}.csv", filename);
    std::ofstream out(out_path);
    if (out) {
      CSVWriter::writeHistogramToStream(out, hist, raw_companion);
    }
  }
}

#ifdef CORRELATION_USE_HDF5
void writeFolderHDF5(const std::filesystem::path &hdf5_dir,
                     const correlation::analysis::DistributionFunctions &dists) {
  std::error_code dir_error;
  std::filesystem::create_directories(hdf5_dir, dir_error);
  const auto out_h5 = hdf5_dir / "analysis.h5";
  HDF5Writer::writeHDF(out_h5.string(), dists);
}
#endif

#ifdef CORRELATION_USE_ARROW
void writeFolderParquet(const std::filesystem::path &arrow_dir,
                        const correlation::analysis::DistributionFunctions &dists,
                        const ExportBundleOptions &options) {
  std::error_code dir_error;
  std::filesystem::create_directories(arrow_dir, dir_error);
  const auto temp_dir = std::filesystem::temp_directory_path() /
                        std::format("correlation_folder_parquet_{}",
                                    std::chrono::steady_clock::now().time_since_epoch().count());
  std::filesystem::create_directories(temp_dir, dir_error);
  const auto base_temp = temp_dir / "parquet_export";
  ArrowWriter::writeAllParquet(base_temp.string(), dists, options.smoothing);

  const auto &all_histograms = dists.getAllHistograms();
  for (const auto &[name, hist] : all_histograms) {
    if (hist.file_suffix.empty() || !isAlgorithmSelected(name, options.selected_algorithms)) {
      continue;
    }
    const auto src_parquet = temp_dir / std::format("parquet_export{}.parquet", hist.file_suffix);
    if (std::filesystem::exists(src_parquet)) {
      std::string filename = hist.file_suffix;
      if (filename.starts_with('_')) {
        filename = filename.substr(1);
      }
      const auto target_parquet = arrow_dir / std::format("{}.parquet", filename);
      std::filesystem::copy_file(src_parquet, target_parquet,
                                 std::filesystem::copy_options::overwrite_existing, dir_error);
    }
  }
  std::filesystem::remove_all(temp_dir, dir_error);
}
#endif

void writeFolderSummary(const std::filesystem::path &dest_dir,
                        const correlation::analysis::DistributionFunctions &dists) {
  const auto summary_path = dest_dir / "summary.txt";
  std::ofstream summary_out(summary_path);
  if (summary_out) {
    summary_out << generateSummaryText(dists);
  }
}

} // namespace

std::expected<void, std::string>
ExportBundleWriter::writeBundle(const std::filesystem::path &zip_path,
                                const correlation::analysis::DistributionFunctions &dists,
                                const ExportBundleOptions &options) {
  // Ensure parent directory exists
  if (zip_path.has_parent_path()) {
    std::error_code dir_error;
    std::filesystem::create_directories(zip_path.parent_path(), dir_error);
    if (dir_error) {
      return std::unexpected(
          std::format("Failed to create destination directory: {}", dir_error.message()));
    }
  }

  utils::ZipArchiveWriter zip;
  auto open_res = zip.open(zip_path);
  if (!open_res) {
    return open_res;
  }

  // 1. Plots
  if (options.include_svg) {
    auto plot_res = packagePlots(zip, dists, options);
    if (!plot_res) {
      return plot_res;
    }
  }

  // 2. CSVs
  if (options.include_csv) {
    auto csv_res = packageCSVs(zip, dists, options);
    if (!csv_res) {
      return csv_res;
    }
  }

  const auto temp_dir = std::filesystem::temp_directory_path() /
                        std::format("correlation_bundle_{}",
                                    std::chrono::steady_clock::now().time_since_epoch().count());
  std::error_code cleanup_error;
  std::filesystem::create_directories(temp_dir, cleanup_error);

  // 3. HDF5
#ifdef CORRELATION_USE_HDF5
  if (options.include_hdf5) {
    auto hdf5_res = packageHDF5(zip, dists, temp_dir);
    if (!hdf5_res) {
      std::filesystem::remove_all(temp_dir, cleanup_error);
      return hdf5_res;
    }
  }
#endif

  // 4. Parquet
#ifdef CORRELATION_USE_ARROW
  if (options.include_parquet) {
    auto parquet_res = packageParquet(zip, dists, options, temp_dir);
    if (!parquet_res) {
      std::filesystem::remove_all(temp_dir, cleanup_error);
      return parquet_res;
    }
  }
#endif

  std::filesystem::remove_all(temp_dir, cleanup_error);

  // 5. Summary
  if (options.include_summary) {
    const std::string summary = generateSummaryText(dists);
    auto summary_res = zip.addFileFromString("summary.txt", summary);
    if (!summary_res) {
      return summary_res;
    }
  }

  return zip.finalize();
}

std::expected<void, std::string>
ExportBundleWriter::writeFolder(const std::filesystem::path &dest_dir,
                                const correlation::analysis::DistributionFunctions &dists,
                                const ExportBundleOptions &options) {
  std::error_code dir_error;
  std::filesystem::create_directories(dest_dir, dir_error);
  if (dir_error) {
    return std::unexpected(
        std::format("Failed to create destination directory: {}", dir_error.message()));
  }

  if (options.include_svg) {
    writeFolderPlots(dest_dir / "plots", dists, options);
  }
  if (options.include_csv) {
    writeFolderCSVs(dest_dir / "csv", dists, options);
  }
#ifdef CORRELATION_USE_HDF5
  if (options.include_hdf5) {
    writeFolderHDF5(dest_dir / "hdf5", dists);
  }
#endif
#ifdef CORRELATION_USE_ARROW
  if (options.include_parquet) {
    writeFolderParquet(dest_dir / "arrow", dists, options);
  }
#endif
  if (options.include_summary) {
    writeFolderSummary(dest_dir, dists);
  }

  return {};
}

} // namespace correlation::writers
