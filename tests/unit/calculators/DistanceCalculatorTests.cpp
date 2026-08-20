#include "analysis/AnalysisTypes.hpp"
#include "calculators/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include <gtest/gtest.h>
#include <vector>

namespace correlation::testing {

using namespace correlation::calculators;
using namespace correlation::core;
using correlation::analysis::BondCutoffMatrix;

TEST(DistanceCalculatorTests, ComputesPairwiseDistancesAndNeighborGraph) {
  // Construct a cubic cell of size 10.0
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});

  // Add 2 atoms: Si (0,0,0) and O (1.5, 0.0, 0.0)
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {1.5, 0.0, 0.0});

  // Setup inputs
  real_t const cutoff_sq = 4.0; // Cutoff distance = 2.0
  BondCutoffMatrix const bond_cutoffs = {
      {{0.36, 4.0}, {0.36, 4.0}}, // Si-Si, Si-O
      {{0.36, 4.0}, {0.36, 4.0}}  // O-Si, O-O
  };

  size_t const num_elements = cell.elements().size();
  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  // Act
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

  // Assert: Check distance histogram bin for 1.5 Å (1.5 / 0.02 = 75)
  // Si (index 0, type 0), O (index 1, type 1)
  ASSERT_EQ(out_histograms.size(), num_elements);
  EXPECT_EQ(out_histograms[0][1][75], 1.0);
  EXPECT_EQ(out_histograms[1][0][75], 1.0);

  // Assert: Check neighbor graph connections
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));

  const auto &neighbors = out_graph.getNeighbors(0);
  ASSERT_EQ(neighbors.size(), 1);
  EXPECT_EQ(neighbors[0].index, 1);
  EXPECT_NEAR(neighbors[0].distance, 1.5, correlation::is_single_precision ? 1e-5 : 1e-9);
}

// --- Extreme / Edge-Case Tests ---

TEST(DistanceCalculatorTests, DistanceAcrossPeriodicBoundary) {
  // Atom near box edge: distance should be computed across PBC
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.5, 5.0, 5.0}); // Near left edge
  cell.addAtom("Si", {9.5, 5.0, 5.0}); // Near right edge
  // PBC distance = 1.0 (not 9.0)

  real_t const cutoff_sq = 4.0; // cutoff = 2.0
  BondCutoffMatrix const bond_cutoffs = {{{0.36, 4.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

  // Distance 1.0 Å in bin 1.0 / 0.02 = 50
  ASSERT_EQ(out_histograms.size(), 1);
  EXPECT_EQ(out_histograms[0][0][50], 1.0);

  // Neighbor graph should reflect the bond
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));
}

TEST(DistanceCalculatorTests, SingleAtomProducesNoDistances) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Ar", {5.0, 5.0, 5.0});

  real_t const cutoff_sq = 25.0;
  BondCutoffMatrix const bond_cutoffs = {{{0.36, 25.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(1);
  DistanceCalculationConfig const config{
      .r_max = 5.0,
      .r_bin_width = 0.02,
      .num_bins = 250,
  };

  // With ignore_periodic_self_interactions = true, a single atom has no pairs
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

  for (real_t const count : out_histograms[0][0]) {
    EXPECT_EQ(count, 0.0);
  }
  EXPECT_TRUE(out_graph.getNeighbors(0).empty());
}

TEST(DistanceCalculatorTests, NonOrthogonalCell) {
  // Triclinic cell
  Cell cell({5.0, 5.0, 5.0, 60.0, 60.0, 60.0});
  cell.addAtom("Ar", {0.0, 0.0, 0.0});
  cell.addAtom("Ar", {2.5, 0.0, 0.0});

  real_t const cutoff_sq = 9.0; // cutoff = 3.0
  BondCutoffMatrix const bond_cutoffs = {{{0.36, 9.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 3.0,
      .r_bin_width = 0.02,
      .num_bins = 150,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

  real_t total_count = 0.0;
  for (real_t const count : out_histograms[0][0]) {
    total_count += count;
  }
  EXPECT_GE(total_count, 1.0);
}

TEST(DistanceCalculatorTests, AtomsOutsideCutoff) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {4.0, 0.0, 0.0}); // Distance = 4.0

  // Cutoff = 2.0, so this pair should NOT be found
  real_t const cutoff_sq = 4.0;
  BondCutoffMatrix const bond_cutoffs = {{{0.36, 4.0}, {0.36, 4.0}}, {{0.36, 4.0}, {0.36, 4.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

  for (real_t const count : out_histograms[0][1]) {
    EXPECT_EQ(count, 0.0);
  }
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
}

TEST(DistanceCalculatorTests, ThrowsOnInvalidCutoff) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  BondCutoffMatrix const bond_cutoffs = {{{0.36, 4.0}}};
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, -1.0, bond_cutoffs, true, out_graph), std::invalid_argument);
  EXPECT_THROW(DistanceCalculator::compute(cell, 0.0, bond_cutoffs, true, out_graph), std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenMinBondCutoffExceedsMax) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // min_sq (4.0) > max_sq (2.25)
  BondCutoffMatrix const invalid_cutoffs = {{{4.0, 2.25}}};
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, invalid_cutoffs, true, out_graph), std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenCutoffBoundsAreNegative) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  NeighborGraph out_graph(1);

  BondCutoffMatrix const neg_min = {{{-0.5, 4.0}}};
  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, neg_min, true, out_graph), std::invalid_argument);

  BondCutoffMatrix const neg_max = {{{0.5, -4.0}}};
  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, neg_max, true, out_graph), std::invalid_argument);
}

TEST(DistanceCalculatorTests, EnforcesMinimumAndMaximumBondCutoffs) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {0.4, 0.0, 0.0}); // Distance = 0.4 Å (below min 0.8 Å)
  cell.addAtom("Si", {1.5, 0.0, 0.0}); // Distance = 1.5 Å (within [0.8, 2.0] Å)
  cell.addAtom("Si", {2.5, 0.0, 0.0}); // Distance = 2.5 Å (above max 2.0 Å, below global cutoff 3.0 Å)

  real_t const cutoff_sq = 9.0; // Global search cutoff = 3.0 Å

  // Bond cutoff window: [0.8 Å, 2.0 Å] -> [0.64 Å², 4.00 Å²]
  BondCutoffMatrix const bond_cutoffs = {{{0.64, 4.00}}};
  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(4);
  DistanceCalculationConfig const config{
      .r_max = 3.0,
      .r_bin_width = 0.02,
      .num_bins = 150,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

  // All 3 pairs are within global cutoff 3.0 Å, so histogram has counts for them
  real_t total_count = 0.0;
  for (real_t const count : out_histograms[0][0]) {
    total_count += count;
  }
  EXPECT_GE(total_count, 3.0);

  // Atom 0 (0.0) -> Atom 1 (0.4 Å): excluded from bonds (0.4 < 0.8 min cutoff)
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));

  // Atom 0 (0.0) -> Atom 2 (1.5 Å): connected as bond (0.8 <= 1.5 <= 2.0)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{2}));

  // Atom 0 (0.0) -> Atom 3 (2.5 Å): excluded from bonds (2.5 > 2.0 max cutoff)
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{3}));
}

} // namespace correlation::testing
