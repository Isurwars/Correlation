/**
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
*/

#include "AppWindow.h"
#include "analysis/PerFrameDistributionStore.hpp"
#include "app/core/AppOptions.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include <gtest/gtest.h>
#include <memory>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

using correlation::real_t;

namespace correlation::app {
namespace {

class Structure3DControllerDistributionTests : public ::testing::Test {
public:
  [[nodiscard]] AppWindow &window() {
    if (!window_.has_value()) {
      throw std::runtime_error("Window not initialized");
    }
    return **window_;
  }

protected:
  void SetUp() override {
#ifndef _WIN32
    setenv("SLINT_BACKEND", "software", 1);
#else
    _putenv_s("SLINT_BACKEND", "software");
#endif
    window_.emplace(AppWindow::create());
  }

private:
  std::optional<slint::ComponentHandle<AppWindow>> window_;
};

TEST_F(Structure3DControllerDistributionTests, SetsStoreAndUpdatesWindowProperties) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto store = std::make_shared<analysis::PerFrameDistributionStore>(3);
  store->setBins(analysis::DistributionMetric::RDF, {1.0, 2.0, 3.0, 4.0});
  store->setMeanCurve(analysis::DistributionMetric::RDF, {0.1, 2.5, 1.2, 0.2});
  store->setFrameCurve(analysis::DistributionMetric::RDF, 0, {0.2, 3.1, 1.5, 0.1});
  store->setFrameCurve(analysis::DistributionMetric::RDF, 1, {0.1, 2.2, 1.1, 0.3});
  store->setFrameCurve(analysis::DistributionMetric::RDF, 2, {0.3, 2.8, 1.3, 0.2});

  controller.setDistributionStore(store);

  EXPECT_TRUE(window().get_structure_3d_dist_has_data());
  EXPECT_TRUE(window().get_structure_3d_dist_hud_visible());
  EXPECT_EQ(window().get_structure_3d_dist_metric(), 0);
  EXPECT_FALSE(std::string(window().get_structure_3d_dist_peak_text().data()).empty());
}

TEST_F(Structure3DControllerDistributionTests, TogglesMetricAndUpdatesWindow) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto store = std::make_shared<analysis::PerFrameDistributionStore>(2);
  store->setBins(analysis::DistributionMetric::PAD, {30.0, 60.0, 90.0, 109.5, 120.0});
  store->setMeanCurve(analysis::DistributionMetric::PAD, {0.0, 0.1, 0.4, 2.8, 0.5});
  store->setFrameCurve(analysis::DistributionMetric::PAD, 0, {0.0, 0.1, 0.5, 3.2, 0.4});

  controller.setDistributionStore(store);

  // Switch to PAD (metric 1)
  controller.handleSelectDistributionMetric(1);

  EXPECT_TRUE(window().get_structure_3d_dist_has_data());
  EXPECT_EQ(window().get_structure_3d_dist_metric(), 1);
  EXPECT_NE(std::string(window().get_structure_3d_dist_peak_text().data()).find("109.5"),
            std::string::npos);
}

TEST_F(Structure3DControllerDistributionTests, TogglesHUDVisibility) {
  const ProgramOptions options;
  Structure3DController controller(window(), options);

  const auto store = std::make_shared<analysis::PerFrameDistributionStore>(1);
  store->setBins(analysis::DistributionMetric::RDF, {1.0, 2.0});
  store->setMeanCurve(analysis::DistributionMetric::RDF, {1.0, 2.0});

  controller.setDistributionStore(store);
  EXPECT_TRUE(window().get_structure_3d_dist_has_data());

  controller.handleToggleDistributionHUD(false);
  EXPECT_FALSE(window().get_structure_3d_dist_has_data());

  controller.handleToggleDistributionHUD(true);
  EXPECT_TRUE(window().get_structure_3d_dist_has_data());
}

} // namespace
} // namespace correlation::app
