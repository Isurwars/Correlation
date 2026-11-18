/**
 * @file AppControllerCalculators.cpp
 * @brief Calculator factory integration, group tree population, and export analysis listing for
 * AppController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/core/AppController.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace correlation::app {

namespace {

bool isDefaultDisabledCalculator(const correlation::calculators::BaseCalculator *calc) noexcept {
  if (calc == nullptr) {
    return false;
  }
  const auto name = calc->getName();
  const auto short_name = calc->getShortName();

  // Enable only Radial (g(r)), Planar Angle Distribution (PAD), Structure Factor (S(Q)), and
  // XRD by default
  const bool is_default_enabled =
      (short_name == "RDF" || name == "Radial Distribution Function" || short_name == "PAD" ||
       name == "Plane Angle Distribution" || short_name == "SQ" ||
       name == "Static Structure Factor" || short_name == "XRD" || name == "X-Ray Diffraction" ||
       name == "X-Ray Diffraction Pattern");
  return !is_default_enabled;
}

using CalculatorPtrList = std::vector<const correlation::calculators::BaseCalculator *>;

/// Calculators bucketed by group, preserving factory registration order of groups.
struct GroupedCalculators {
  std::vector<std::string> order;
  std::map<std::string, CalculatorPtrList, std::less<>> members;
};

[[nodiscard]] GroupedCalculators groupRegisteredCalculators() {
  GroupedCalculators grouped;
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  for (const auto &calc : calculators) {
    const std::string grp(calc->getGroup());
    if (!grouped.members.contains(grp)) {
      grouped.order.push_back(grp);
    }
    grouped.members[grp].push_back(calc.get());
  }
  return grouped;
}

[[nodiscard]] bool isVdosEnabled(const ProgramOptions &opts) noexcept {
  const auto vdos_it = opts.active_calculators.find("Vibrational Density of States");
  if (vdos_it != opts.active_calculators.end() && vdos_it->second) {
    return true;
  }
  const auto vdos_short_it = opts.active_calculators.find("VDOS");
  if (vdos_short_it != opts.active_calculators.end() && vdos_short_it->second) {
    return true;
  }
  const auto vdos_alias_it = opts.active_calculators.find("vDoS");
  return vdos_alias_it != opts.active_calculators.end() && vdos_alias_it->second;
}

[[nodiscard]] bool resolveCalculatorEnabled(const correlation::calculators::BaseCalculator &calc,
                                            const ProgramOptions &opts) {
  if (const auto calc_it = opts.active_calculators.find(std::string(calc.getName()));
      calc_it != opts.active_calculators.end()) {
    return calc_it->second;
  }
  if (const auto short_it = opts.active_calculators.find(std::string(calc.getShortName()));
      short_it != opts.active_calculators.end()) {
    return short_it->second;
  }
  if ((calc.getShortName() == "VACF" || calc.getName() == "Velocity Autocorrelation Function") &&
      isVdosEnabled(opts)) {
    return true;
  }
  return !isDefaultDisabledCalculator(&calc);
}

[[nodiscard]] CalculatorGroup buildCalculatorGroup(const std::string &group_name,
                                                   const CalculatorPtrList &calcs,
                                                   const ProgramOptions &opts, bool expanded) {
  const auto calcs_model = std::make_shared<slint::VectorModel<CalculatorInfo>>();
  int enabled_count = 0;
  for (const auto *calc : calcs) {
    const bool enabled = resolveCalculatorEnabled(*calc, opts);
    enabled_count += enabled ? 1 : 0;
    CalculatorInfo info;
    info.id = slint::SharedString(calc->getName());
    info.name = slint::SharedString(calc->getName());
    info.short_name = slint::SharedString(calc->getShortName());
    info.description = slint::SharedString(calc->getDescription());
    info.enabled = enabled;
    calcs_model->push_back(info);
  }
  CalculatorGroup group;
  group.name = slint::SharedString(group_name);
  group.expanded = expanded;
  group.enabled_count = enabled_count;
  group.total_count = static_cast<int>(calcs.size());
  group.calculators = calcs_model;
  return group;
}

/// True when @p model already holds exactly the groups in @p order (same names, same order).
[[nodiscard]] bool hasSameGroupLayout(const std::shared_ptr<slint::Model<CalculatorGroup>> &model,
                                      const std::vector<std::string> &order) {
  if (model == nullptr || model->row_count() != order.size()) {
    return false;
  }
  for (size_t idx = 0; idx < order.size(); ++idx) {
    const auto row = model->row_data(idx);
    if (!row.has_value() || std::string_view(row->name.data()) != order[idx]) {
      return false;
    }
  }
  return true;
}

} // namespace

bool AppController::isCalculatorGroupExpanded(std::string_view group_name,
                                              size_t group_index) const {
  if (const auto found_it = calculator_group_expanded_.find(group_name);
      found_it != calculator_group_expanded_.end()) {
    return found_it->second;
  }
  return group_name == "Radial" || group_index == 0;
}

void AppController::populateCalculatorGroups() {
  const auto grouped = groupRegisteredCalculators();

  std::vector<CalculatorGroup> groups;
  groups.reserve(grouped.order.size());
  int total_count = 0;
  int total_enabled = 0;
  for (size_t idx = 0; idx < grouped.order.size(); ++idx) {
    const auto &grp_name = grouped.order[idx];
    auto group = buildCalculatorGroup(grp_name, grouped.members.at(grp_name), options_,
                                      isCalculatorGroupExpanded(grp_name, idx));
    total_count += group.total_count;
    total_enabled += group.enabled_count;
    groups.push_back(std::move(group));
  }

  // Update rows in place when the layout is unchanged so that Slint keeps the existing
  // delegates (open sub-trees, scroll offset) instead of rebuilding the whole tree.
  const auto existing = window_.get_calculator_groups();
  if (hasSameGroupLayout(existing, grouped.order)) {
    for (size_t idx = 0; idx < groups.size(); ++idx) {
      existing->set_row_data(idx, groups[idx]);
    }
  } else {
    const auto groups_model = std::make_shared<slint::VectorModel<CalculatorGroup>>();
    for (const auto &group : groups) {
      groups_model->push_back(group);
    }
    window_.set_calculator_groups(groups_model);
  }
  window_.set_total_calculator_count(total_count);
  window_.set_enabled_calculator_count(total_enabled);
}

void AppController::setCalculatorsEnabled(bool enabled, std::string_view group_filter) {
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  if (group_filter.empty() && !enabled) {
    options_.active_calculators.clear();
    for (const auto &calc : calculators) {
      options_.active_calculators[std::string(calc->getName())] = false;
      if (!calc->getShortName().empty()) {
        options_.active_calculators[std::string(calc->getShortName())] = false;
      }
    }
  } else {
    for (const auto &calc : calculators) {
      if (!group_filter.empty() && calc->getGroup() != group_filter) {
        continue;
      }
      options_.active_calculators[std::string(calc->getName())] = enabled;
      if (!calc->getShortName().empty()) {
        options_.active_calculators[std::string(calc->getShortName())] = enabled;
      }
    }
    if (isVdosEnabled(options_)) {
      options_.active_calculators["Velocity Autocorrelation Function"] = true;
      options_.active_calculators["VACF"] = true;
    }
  }
  populateCalculatorGroups();
  updateActiveGroupFlags();
}

void AppController::setCalculatorGroupsExpanded(bool expanded, std::string_view group_filter) {
  if (!group_filter.empty()) {
    calculator_group_expanded_.insert_or_assign(std::string(group_filter), expanded);
  } else {
    for (const auto &grp_name : groupRegisteredCalculators().order) {
      calculator_group_expanded_.insert_or_assign(grp_name, expanded);
    }
  }
  populateCalculatorGroups();
}

void AppController::populateExportAnalyses() {
  const auto names = dispatcher_.getAvailableHistogramNames();
  const auto model = std::make_shared<slint::VectorModel<ExportAnalysisItem>>();
  for (const auto &name : names) {
    if (name.ends_with("_raw")) {
      continue;
    }
    const auto *hist = dispatcher_.getHistogram(name);
    const std::string label = (hist != nullptr && !hist->title.empty()) ? hist->title : name;
    ExportAnalysisItem item;
    item.id = slint::SharedString(name);
    item.name = slint::SharedString(label);
    item.enabled = true;
    model->push_back(item);
  }
  window_.set_export_analysis_items(model);
  window_.set_export_all_analyses(true);
}

} // namespace correlation::app
