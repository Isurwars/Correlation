/**
 * @file PlotRenderPipelineTests.cpp
 * @brief Unit tests for PlotRenderPipeline.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotRenderPipeline.hpp"
#include <condition_variable>
#include <gtest/gtest.h>
#include <mutex>

namespace correlation::app {
namespace {

using namespace std::chrono_literals;

TEST(PlotRenderPipelineTests, InitialCacheStateReturnsMiss) {
  const PlotRenderPipeline pipeline([](const slint::Image &) {}, false);
  const correlation::plotters::PlotConfig config;
  const correlation::plotters::HoverInfo hover;

  EXPECT_FALSE(pipeline.isCacheHit(0, config, hover, 0));
}

TEST(PlotRenderPipelineTests, CacheHitDetectionMatchesIdenticalTask) {
  PlotRenderPipeline pipeline([](const slint::Image &) {}, false);
  correlation::plotters::PlotConfig config;
  config.width = 800.0;
  config.height = 600.0;
  config.theme = correlation::plotters::PlotConfig::Theme::Dark;

  correlation::plotters::HoverInfo hover;
  hover.active = false;

  RenderTask task;
  task.index = 2;
  task.config = config;
  task.hover = hover;
  task.pinned_runs_count = 3;

  pipeline.submitTask(task);

  EXPECT_TRUE(pipeline.isCacheHit(2, config, hover, 3));
  EXPECT_FALSE(pipeline.isCacheHit(1, config, hover, 3)); // Different index
  EXPECT_FALSE(pipeline.isCacheHit(2, config, hover, 2)); // Different pinned count

  auto modified_config = config;
  modified_config.theme = correlation::plotters::PlotConfig::Theme::Light;
  EXPECT_FALSE(pipeline.isCacheHit(2, modified_config, hover, 3)); // Different theme
}

TEST(PlotRenderPipelineTests, InvalidateCacheForcesMiss) {
  PlotRenderPipeline pipeline([](const slint::Image &) {}, false);
  const correlation::plotters::PlotConfig config;
  const correlation::plotters::HoverInfo hover;

  RenderTask task;
  task.index = 1;
  task.config = config;
  task.hover = hover;
  task.pinned_runs_count = 0;

  pipeline.submitTask(task);
  EXPECT_TRUE(pipeline.isCacheHit(1, config, hover, 0));

  pipeline.invalidateCache();
  EXPECT_FALSE(pipeline.isCacheHit(1, config, hover, 0));
}

TEST(PlotRenderPipelineTests, TaskExecutionRendersImageAndInvokesCallback) {
  std::mutex cv_mutex;
  std::condition_variable cond_var;
  bool callback_invoked{false};
  slint::Image rendered_image;

  PlotRenderPipeline pipeline(
      [&](const slint::Image &img) {
        {
          const std::scoped_lock lock(cv_mutex);
          callback_invoked = true;
          rendered_image = img;
        }
        cond_var.notify_one();
      },
      false);

  RenderTask task;
  task.index = 0;
  task.active_hist.bins = {1.0, 2.0, 3.0};
  task.active_hist.partials["Total"] = {10.0, 20.0, 5.0};
  task.config.width = 400.0;
  task.config.height = 300.0;

  pipeline.submitTask(std::move(task));

  std::unique_lock lock(cv_mutex);
  const bool received = cond_var.wait_for(lock, 2s, [&] { return callback_invoked; });
  EXPECT_TRUE(received);
  EXPECT_GT(rendered_image.size().width, 0U);
  EXPECT_GT(rendered_image.size().height, 0U);
}

TEST(PlotRenderPipelineTests, StopShutsDownCleanly) {
  PlotRenderPipeline pipeline([](const slint::Image &) {}, false);
  pipeline.stop();
  SUCCEED();
}

} // namespace
} // namespace correlation::app
