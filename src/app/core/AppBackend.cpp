/**
 * @file AppBackend.cpp
 * @brief Implementation of AppBackend facade delegating to specialized services.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/core/AppBackend.hpp"
#include <iostream>

namespace correlation::app {

AppBackend::AppBackend() = default;

std::string AppBackend::loadFile(const std::string &path) {
  auto res = loader_.loadFile(path, progress_callback_);
  if (!res) {
    throw std::runtime_error(res.error());
  }
  options_.input_file = path;
  options_.output_file_base = path;
  return res.value();
}

correlation::analysis::BondCutoffMatrix AppBackend::getRecommendedBondCutoffs() {
  return cutoff_service_.getRecommendedBondCutoffs(loader_.trajectoryMut());
}

void AppBackend::setBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs) {
  cutoff_service_.setBondCutoffs(loader_.trajectoryMut(), cutoffs);
  options_.bond_cutoffs = cutoffs;
}

correlation::analysis::BondCutoffMatrix AppBackend::applyScaledBondCutoffs(real_t scale_factor) {
  auto cutoffs = cutoff_service_.applyScaledBondCutoffs(loader_.trajectoryMut(), scale_factor);
  options_.bond_cutoffs = cutoffs;
  return cutoffs;
}

correlation::analysis::BondCutoffMatrix AppBackend::setUniformBondCutoff(real_t min_cutoff,
                                                                         real_t max_cutoff) {
  auto cutoffs = cutoff_service_.setUniformBondCutoff(loader_.cell(), loader_.trajectoryMut(),
                                                      min_cutoff, max_cutoff);
  options_.bond_cutoffs = cutoffs;
  return cutoffs;
}

std::expected<void, std::string> AppBackend::runAnalysis() {
  if (loader_.trajectory() == nullptr || loader_.getFrameCount() == 0) {
    std::string const err = AppDefaults::MSG_ANALYSIS_ABORTED;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }
  dispatcher_.setProgressCallback(progress_callback_);
  return dispatcher_.runAnalysis(*loader_.trajectoryMut(), options_);
}

std::expected<void, std::string> AppBackend::writeFiles() {
  return dispatcher_.writeFiles(options_);
}

} // namespace correlation::app
