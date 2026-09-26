// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/formatters/BondCutoffMapper.hpp"
#include <gtest/gtest.h>

namespace {

using correlation::app::BondCutoffMapper;
using correlation::app::CutoffEntry;
using correlation::core::Element;

TEST(BondCutoffMapperTests, DefaultCutoffEntriesGeneration) {
  std::vector<Element> elements = {
      Element{.symbol = "Si", .id = correlation::core::ElementID{0}},
      Element{.symbol = "O", .id = correlation::core::ElementID{1}},
  };

  const auto entries = BondCutoffMapper::createDefaultCutoffEntries(elements);
  // Pair combinations for 2 elements: (Si-Si, Si-O, O-O) = 3
  ASSERT_EQ(entries.size(), 3U);

  EXPECT_EQ(entries[0].element1, "Si");
  EXPECT_EQ(entries[0].element2, "Si");
  EXPECT_FALSE(entries[0].min_distance.empty());
  EXPECT_FALSE(entries[0].max_distance.empty());

  EXPECT_EQ(entries[1].element1, "Si");
  EXPECT_EQ(entries[1].element2, "O");

  EXPECT_EQ(entries[2].element1, "O");
  EXPECT_EQ(entries[2].element2, "O");
}

TEST(BondCutoffMapperTests, ParseCutoffMatrixSymmetryAndValues) {
  std::vector<Element> elements = {
      Element{.symbol = "Si", .id = correlation::core::ElementID{0}},
      Element{.symbol = "O", .id = correlation::core::ElementID{1}},
  };

  std::vector<CutoffEntry> entries = {
      CutoffEntry{.element1 = "Si", .element2 = "Si", .min_distance = "1.0", .max_distance = "3.0"},
      CutoffEntry{.element1 = "Si", .element2 = "O", .min_distance = "1.2", .max_distance = "2.2"},
      CutoffEntry{.element1 = "O", .element2 = "O", .min_distance = "1.1", .max_distance = "2.8"},
  };

  const auto matrix = BondCutoffMapper::parseCutoffMatrix(entries, elements);
  ASSERT_EQ(matrix.size(), 2U);
  ASSERT_EQ(matrix[0].size(), 2U);
  ASSERT_EQ(matrix[1].size(), 2U);

  // Si-Si (idx 0, 0)
  EXPECT_DOUBLE_EQ(matrix[0][0].min_sq, 1.0);
  EXPECT_DOUBLE_EQ(matrix[0][0].max_sq, 9.0);

  // Si-O (idx 0, 1) and O-Si (idx 1, 0) symmetric
  EXPECT_NEAR(matrix[0][1].min_sq, 1.44, 1e-4);
  EXPECT_NEAR(matrix[0][1].max_sq, 4.84, 1e-4);
  EXPECT_DOUBLE_EQ(matrix[0][1].min_sq, matrix[1][0].min_sq);
  EXPECT_DOUBLE_EQ(matrix[0][1].max_sq, matrix[1][0].max_sq);

  // O-O (idx 1, 1)
  EXPECT_NEAR(matrix[1][1].min_sq, 1.21, 1e-4);
  EXPECT_NEAR(matrix[1][1].max_sq, 7.84, 1e-4);
}

TEST(BondCutoffMapperTests, ParseCutoffMatrixInvalidDistanceHandledGracefully) {
  std::vector<Element> elements = {
      Element{.symbol = "C", .id = correlation::core::ElementID{0}},
  };

  std::vector<CutoffEntry> entries = {
      CutoffEntry{.element1 = "C",
                  .element2 = "C",
                  .min_distance = "invalid_min",
                  .max_distance = "not_a_number"},
  };

  const auto matrix = BondCutoffMapper::parseCutoffMatrix(entries, elements);
  ASSERT_EQ(matrix.size(), 1U);
  EXPECT_DOUBLE_EQ(matrix[0][0].min_sq, 0.0);
  EXPECT_DOUBLE_EQ(matrix[0][0].max_sq, 0.0);
}

} // namespace
