/**
 * @file BondCutoffController.cpp
 * @brief Implementation of BondCutoffController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/BondCutoffController.hpp"
#include "app/AppBackend.hpp"
#include "app/BondCutoffMapper.hpp"
#include "physics/PhysicalData.hpp"

#include <cmath>
#include <format>
#include <memory>
#include <string>
#include <vector>

namespace correlation::app {

BondCutoffController::BondCutoffController(AppWindow &window, AppBackend &backend)
    : window_(window), backend_(backend) {}

void BondCutoffController::setBondCutoffs() {
  if (backend_.cell() == nullptr) {
    return;
  }

  const auto &elements = backend_.cell()->elements();
  const auto entries = BondCutoffMapper::createDefaultCutoffEntries(elements);

  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (const auto &entry : entries) {
    slint_cutoffs->push_back({
        .element1 = slint::SharedString(entry.element1),
        .element2 = slint::SharedString(entry.element2),
        .min_distance = slint::SharedString(entry.min_distance),
        .max_distance = slint::SharedString(entry.max_distance),
    });
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  backend_.setBondCutoffs(getBondCutoffs());
}

correlation::analysis::BondCutoffMatrix BondCutoffController::getBondCutoffs() const {
  auto slint_cutoffs = window_.get_bond_cutoffs();
  if (backend_.cell() == nullptr || slint_cutoffs == nullptr) {
    return {};
  }

  const auto &elements = backend_.cell()->elements();
  std::vector<CutoffEntry> entries;
  entries.reserve(slint_cutoffs->row_count());

  for (size_t k = 0; k < slint_cutoffs->row_count(); ++k) {
    auto maybe_item = slint_cutoffs->row_data(k);
    if (!maybe_item.has_value()) {
      continue;
    }
    const auto &item = maybe_item.value();
    entries.push_back(CutoffEntry{
        .element1 = std::string(item.element1.data()),
        .element2 = std::string(item.element2.data()),
        .min_distance = std::string(item.min_distance.data()),
        .max_distance = std::string(item.max_distance.data()),
    });
  }

  return BondCutoffMapper::parseCutoffMatrix(entries, elements);
}

void BondCutoffController::applyScaledCutoffs(float scale_factor) {
  if (backend_.cell() == nullptr || scale_factor <= 0.0F) {
    return;
  }
  const auto scaled_cutoffs = backend_.applyScaledBondCutoffs(static_cast<real_t>(scale_factor));
  const auto &elements = backend_.cell()->elements();
  const auto num_elements = elements.size();

  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t min_d = std::sqrt(scaled_cutoffs[i][j].min_sq);
      const real_t max_d = std::sqrt(scaled_cutoffs[i][j].max_sq);
      slint_cutoffs->push_back(BondCutoff{
          .element1 = slint::SharedString(elements[i].symbol),
          .element2 = slint::SharedString(elements[j].symbol),
          .min_distance = slint::SharedString(std::format("{:.2f}", min_d)),
          .max_distance = slint::SharedString(std::format("{:.2f}", max_d)),
      });
    }
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  backend_.setBondCutoffs(getBondCutoffs());
}

void BondCutoffController::setUniformCutoff(float max_cutoff) {
  if (backend_.cell() == nullptr || max_cutoff <= 0.0F) {
    return;
  }
  const auto uniform_cutoffs = backend_.setUniformBondCutoff(0.0, static_cast<real_t>(max_cutoff));
  const auto &elements = backend_.cell()->elements();
  const auto num_elements = elements.size();

  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t min_d = std::sqrt(uniform_cutoffs[i][j].min_sq);
      const real_t max_d = std::sqrt(uniform_cutoffs[i][j].max_sq);
      slint_cutoffs->push_back(BondCutoff{
          .element1 = slint::SharedString(elements[i].symbol),
          .element2 = slint::SharedString(elements[j].symbol),
          .min_distance = slint::SharedString(std::format("{:.2f}", min_d)),
          .max_distance = slint::SharedString(std::format("{:.2f}", max_d)),
      });
    }
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  backend_.setBondCutoffs(getBondCutoffs());
}

namespace {

[[nodiscard]] real_t safeGetCovalentRadius(const std::string &symbol) {
  try {
    return physics::getCovalentRadius(symbol);
  } catch (const std::out_of_range &) {
    return static_cast<real_t>(1.5);
  }
}

struct CovalentCutoffInput {
  std::string symbol_a;
  std::string symbol_b;
  real_t radius_a = 0.0;
  real_t radius_b = 0.0;
  size_t row_idx = 0;
  float factor = 1.0F;
  FactorBound bound = FactorBound::Min;
};

[[nodiscard]] BondCutoff
createCovalentCutoffRow(const CovalentCutoffInput &input,
                        const std::shared_ptr<slint::Model<BondCutoff>> &slint_cutoffs) {
  const real_t sum_radii = input.radius_a + input.radius_b;
  slint::SharedString min_str;
  slint::SharedString max_str;

  if (input.bound == FactorBound::Min) {
    const real_t min_dist = sum_radii * static_cast<real_t>(input.factor);
    min_str = slint::SharedString(std::format("{:.2f}", min_dist));
    if (slint_cutoffs != nullptr && input.row_idx < slint_cutoffs->row_count()) {
      const auto opt = slint_cutoffs->row_data(input.row_idx);
      max_str = opt.has_value() ? opt->max_distance
                                : slint::SharedString(std::format(
                                      "{:.2f}", sum_radii * AppDefaults::BOND_MAX_FACTOR));
    } else {
      max_str =
          slint::SharedString(std::format("{:.2f}", sum_radii * AppDefaults::BOND_MAX_FACTOR));
    }
  } else {
    const real_t max_dist = sum_radii * static_cast<real_t>(input.factor);
    max_str = slint::SharedString(std::format("{:.2f}", max_dist));
    if (slint_cutoffs != nullptr && input.row_idx < slint_cutoffs->row_count()) {
      const auto opt = slint_cutoffs->row_data(input.row_idx);
      min_str = opt.has_value() ? opt->min_distance
                                : slint::SharedString(std::format(
                                      "{:.2f}", sum_radii * AppDefaults::BOND_MIN_FACTOR));
    } else {
      min_str =
          slint::SharedString(std::format("{:.2f}", sum_radii * AppDefaults::BOND_MIN_FACTOR));
    }
  }

  return BondCutoff{
      .element1 = slint::SharedString(input.symbol_a),
      .element2 = slint::SharedString(input.symbol_b),
      .min_distance = min_str,
      .max_distance = max_str,
  };
}

} // namespace

void BondCutoffController::applyCovalentFactor(float factor, FactorBound bound) {
  if (backend_.cell() == nullptr || factor <= 0.0F) {
    return;
  }
  const auto &elements = backend_.cell()->elements();
  const auto num_elements = elements.size();
  const auto slint_cutoffs = window_.get_bond_cutoffs();

  auto new_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  size_t row_index = 0;
  for (size_t i = 0; i < num_elements; ++i) {
    const real_t radius_a = safeGetCovalentRadius(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const real_t radius_b = safeGetCovalentRadius(elements[j].symbol);
      new_cutoffs->push_back(createCovalentCutoffRow(
          CovalentCutoffInput{
              .symbol_a = elements[i].symbol,
              .symbol_b = elements[j].symbol,
              .radius_a = radius_a,
              .radius_b = radius_b,
              .row_idx = row_index,
              .factor = factor,
              .bound = bound,
          },
          slint_cutoffs));
      ++row_index;
    }
  }

  window_.set_bond_cutoffs(new_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  backend_.setBondCutoffs(getBondCutoffs());
}

void BondCutoffController::applyGlobalCutoff(float global_cutoff) {
  setUniformCutoff(global_cutoff);
}

} // namespace correlation::app
