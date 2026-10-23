/**
 * @file SteinhardtLRUCacheTests.cpp
 * @brief Unit tests for SteinhardtLRUCache and PhaseSearchTypes.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PhaseSearchTypes.hpp"
#include "app/viewmodel/SteinhardtLRUCache.hpp"

#include <gtest/gtest.h>
#include <thread>
#include <vector>

namespace correlation::app {
namespace {

using SteinhardtParams = calculators::SteinhardtCalculator::SteinhardtParams;

std::shared_ptr<SteinhardtParams> makeDummyParams(std::size_t num_atoms, real_t base_q6) {
  auto params = std::make_shared<SteinhardtParams>();
  params->Q6_bar.assign(num_atoms, base_q6);
  params->Q4_bar.assign(num_atoms, base_q6 * static_cast<real_t>(0.5));
  params->Q6.assign(num_atoms, base_q6);
  params->Q4.assign(num_atoms, base_q6 * static_cast<real_t>(0.5));
  params->W6_hat.assign(num_atoms, -0.01);
  params->W4_hat.assign(num_atoms, -0.15);
  return params;
}

TEST(SteinhardtLRUCacheTests, InitialStateAndCapacity) {
  const SteinhardtLRUCache cache(16);
  EXPECT_EQ(cache.capacity(), 16UL);
  EXPECT_EQ(cache.size(), 0UL);
  EXPECT_FALSE(cache.contains(0));
  EXPECT_EQ(cache.peek(0), nullptr);
}

TEST(SteinhardtLRUCacheTests, PutAndRetrieve) {
  SteinhardtLRUCache cache(4);
  const auto params = makeDummyParams(10, 0.575);

  cache.put(0, params);
  EXPECT_EQ(cache.size(), 1UL);
  EXPECT_TRUE(cache.contains(0));

  const auto retrieved = cache.get(0);
  ASSERT_NE(retrieved, nullptr);
  EXPECT_EQ(retrieved->Q6_bar.size(), 10UL);
  EXPECT_NEAR(retrieved->Q6_bar[0], 0.575, 1e-5);
}

TEST(SteinhardtLRUCacheTests, LRUEvictionOrder) {
  // Cache with capacity 3
  SteinhardtLRUCache cache(3);
  cache.put(1, makeDummyParams(1, 0.1));
  cache.put(2, makeDummyParams(1, 0.2));
  cache.put(3, makeDummyParams(1, 0.3));

  EXPECT_EQ(cache.size(), 3UL);

  // Access frame 1 (MRU: 1, then 3, then 2 oldest)
  EXPECT_NE(cache.get(1), nullptr);

  // Insert frame 4 -> frame 2 should be evicted
  cache.put(4, makeDummyParams(1, 0.4));

  EXPECT_EQ(cache.size(), 3UL);
  EXPECT_TRUE(cache.contains(1));
  EXPECT_FALSE(cache.contains(2)); // Evicted!
  EXPECT_TRUE(cache.contains(3));
  EXPECT_TRUE(cache.contains(4));
}

TEST(SteinhardtLRUCacheTests, PeekDoesNotPromoteToMRU) {
  SteinhardtLRUCache cache(2);
  cache.put(10, makeDummyParams(1, 0.1));
  cache.put(20, makeDummyParams(1, 0.2));

  // Peek frame 10 (should NOT make it MRU; frame 10 remains oldest)
  EXPECT_NE(cache.peek(10), nullptr);

  // Put frame 30 -> frame 10 should be evicted because peek didn't update order
  cache.put(30, makeDummyParams(1, 0.3));

  EXPECT_FALSE(cache.contains(10));
  EXPECT_TRUE(cache.contains(20));
  EXPECT_TRUE(cache.contains(30));
}

TEST(SteinhardtLRUCacheTests, ClearPurgesAllEntries) {
  SteinhardtLRUCache cache(4);
  cache.put(1, makeDummyParams(5, 0.5));
  cache.put(2, makeDummyParams(5, 0.6));
  EXPECT_EQ(cache.size(), 2UL);

  cache.clear();
  EXPECT_EQ(cache.size(), 0UL);
  EXPECT_FALSE(cache.contains(1));
  EXPECT_FALSE(cache.contains(2));
}

TEST(SteinhardtLRUCacheTests, ConcurrentAccessThreadSafety) {
  SteinhardtLRUCache cache(8);
  std::vector<std::jthread> threads;
  threads.reserve(8);

  for (std::size_t thread_id = 0; thread_id < 8; ++thread_id) {
    threads.emplace_back([&cache, thread_id] {
      for (std::size_t iter = 0; iter < 100; ++iter) {
        const auto frame = (thread_id * 10 + iter) % 16;
        cache.put(frame, makeDummyParams(4, 0.5));
        static_cast<void>(cache.contains(frame));
        static_cast<void>(cache.peek(frame));
        static_cast<void>(cache.get(frame));
      }
    });
  }

  for (auto &thread_instance : threads) {
    if (thread_instance.joinable()) {
      thread_instance.join();
    }
  }

  EXPECT_LE(cache.size(), cache.capacity());
}

TEST(PhaseSearchTypesTests, CanonicalPresetConfigurations) {
  const auto liquid_cfg = getPresetConfig(CrystalPreset::Liquid);
  EXPECT_EQ(liquid_cfg.descriptor, PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(liquid_cfg.min_threshold, 0.10F);
  EXPECT_FLOAT_EQ(liquid_cfg.max_threshold, 0.28F);

  const auto fcc_cfg = getPresetConfig(CrystalPreset::FCC);
  EXPECT_EQ(fcc_cfg.descriptor, PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(fcc_cfg.min_threshold, 0.46F);
  EXPECT_FLOAT_EQ(fcc_cfg.max_threshold, 0.58F);

  const auto hcp_cfg = getPresetConfig(CrystalPreset::HCP);
  EXPECT_EQ(hcp_cfg.descriptor, PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(hcp_cfg.min_threshold, 0.44F);
  EXPECT_FLOAT_EQ(hcp_cfg.max_threshold, 0.52F);

  const auto bcc_cfg = getPresetConfig(CrystalPreset::BCC);
  EXPECT_EQ(bcc_cfg.descriptor, PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(bcc_cfg.min_threshold, 0.35F);
  EXPECT_FLOAT_EQ(bcc_cfg.max_threshold, 0.44F);

  const auto ico_cfg = getPresetConfig(CrystalPreset::Icosahedral);
  EXPECT_EQ(ico_cfg.descriptor, PhaseDescriptor::W6Hat);
  EXPECT_FLOAT_EQ(ico_cfg.min_threshold, -0.17F);
  EXPECT_FLOAT_EQ(ico_cfg.max_threshold, -0.14F);

  const auto custom_cfg = getPresetConfig(CrystalPreset::Custom);
  EXPECT_EQ(custom_cfg.descriptor, PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(custom_cfg.min_threshold, 0.0F);
  EXPECT_FLOAT_EQ(custom_cfg.max_threshold, 1.0F);
}

TEST(PhaseSearchTypesTests, StringConversionsAndDomains) {
  EXPECT_EQ(descriptorToString(PhaseDescriptor::Q6Bar), "q̄₆");
  EXPECT_EQ(descriptorToString(PhaseDescriptor::W6Hat), "ŵ₆");
  EXPECT_EQ(descriptorToString(PhaseDescriptor::W4Hat), "ŵ₄");
  EXPECT_EQ(descriptorToString(PhaseDescriptor::LocalEntropy), "s₂ (LEF)");
  EXPECT_EQ(descriptorToString(PhaseDescriptor::TWF_Nconn), "N_conn (TWF)");

  EXPECT_EQ(presetToString(CrystalPreset::FCC), "FCC");
  EXPECT_EQ(presetToString(CrystalPreset::Icosahedral), "Icosahedral");

  const auto [q_min, q_max] = getDescriptorDomain(PhaseDescriptor::Q6Bar);
  EXPECT_FLOAT_EQ(q_min, 0.0F);
  EXPECT_FLOAT_EQ(q_max, 1.0F);

  const auto [w_min, w_max] = getDescriptorDomain(PhaseDescriptor::W6Hat);
  EXPECT_FLOAT_EQ(w_min, -0.2F);
  EXPECT_FLOAT_EQ(w_max, 0.2F);

  const auto [w4_min, w4_max] = getDescriptorDomain(PhaseDescriptor::W4Hat);
  EXPECT_FLOAT_EQ(w4_min, -0.2F);
  EXPECT_FLOAT_EQ(w4_max, 0.2F);

  const auto [s2_min, s2_max] = getDescriptorDomain(PhaseDescriptor::LocalEntropy);
  EXPECT_FLOAT_EQ(s2_min, -15.0F);
  EXPECT_FLOAT_EQ(s2_max, 0.0F);

  const auto [nconn_min, nconn_max] = getDescriptorDomain(PhaseDescriptor::TWF_Nconn);
  EXPECT_FLOAT_EQ(nconn_min, 0.0F);
  EXPECT_FLOAT_EQ(nconn_max, 16.0F);
}

TEST(PhaseDescriptorsLRUCacheTests, StoreAndRetrieveUnifiedDescriptors) {
  PhaseDescriptorsLRUCache cache(4);
  auto desc = std::make_shared<FramePhaseDescriptors>();
  desc->steinhardt = makeDummyParams(10, 0.55);
  desc->local_entropy = {-6.15, -6.10, -6.20};
  desc->twf_nconn = {12.0, 11.0, 12.0};

  cache.put(0, desc);
  EXPECT_TRUE(cache.contains(0));
  const auto retrieved = cache.get(0);
  ASSERT_NE(retrieved, nullptr);
  ASSERT_NE(retrieved->steinhardt, nullptr);
  EXPECT_EQ(retrieved->steinhardt->Q6_bar.size(), 10UL);
  EXPECT_EQ(retrieved->local_entropy.size(), 3UL);
  EXPECT_NEAR(retrieved->local_entropy[0], -6.15, 1e-5);
  ASSERT_EQ(retrieved->twf_nconn.size(), 3UL);
  EXPECT_FLOAT_EQ(static_cast<float>(retrieved->twf_nconn[0]), 12.0F);
}

} // namespace
} // namespace correlation::app
