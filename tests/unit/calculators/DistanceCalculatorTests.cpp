#include "analysis/AnalysisTypes.hpp"
#include "calculators/spatial/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include "../../CrystalTestHelper.hpp"

#include <algorithm>
#include <gtest/gtest.h>
#include <numbers>
#include <random>
#include <tbb/global_control.h>
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
      {{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}}, // Si-Si, Si-O
      {{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}}, // O-Si, O-O
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
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 4.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 25.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(1);
  DistanceCalculationConfig const config{
      .r_max = 5.0,
      .r_bin_width = 0.02,
      .num_bins = 250,
  };

  // With ignore_periodic_self_interactions = true, a single atom has no pairs
  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 9.0}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 3.0,
      .r_bin_width = 0.02,
      .num_bins = 150,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

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
  BondCutoffMatrix const bond_cutoffs = {
      {{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}},
      {{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}},
  };

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(2);
  DistanceCalculationConfig const config{
      .r_max = 2.0,
      .r_bin_width = 0.02,
      .num_bins = 100,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

  for (real_t const count : out_histograms[0][1]) {
    EXPECT_EQ(count, 0.0);
  }
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
}

TEST(DistanceCalculatorTests, ThrowsOnInvalidCutoff) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 4.0}}};
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, -1.0, bond_cutoffs, true, out_graph),
               std::invalid_argument);
  EXPECT_THROW(DistanceCalculator::compute(cell, 0.0, bond_cutoffs, true, out_graph),
               std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenMinBondCutoffExceedsMax) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // min_sq (4.0) > max_sq (2.25)
  BondCutoffMatrix const invalid_cutoffs = {{{.min_sq = 4.0, .max_sq = 2.25}}};
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, invalid_cutoffs, true, out_graph),
               std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenCutoffBoundsAreNegative) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  NeighborGraph out_graph(1);

  BondCutoffMatrix const neg_min = {{{.min_sq = -0.5, .max_sq = 4.0}}};
  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, neg_min, true, out_graph),
               std::invalid_argument);

  BondCutoffMatrix const neg_max = {{{.min_sq = 0.5, .max_sq = -4.0}}};
  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, neg_max, true, out_graph),
               std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenBondCutoffExceedsSearchCutoff) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  NeighborGraph out_graph(1);

  // max_sq (16.0) > cutoff_sq (9.0)
  BondCutoffMatrix const exceeds_cutoff = {{{.min_sq = 0.5, .max_sq = 16.0}}};
  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, exceeds_cutoff, true, out_graph),
               std::invalid_argument);
}

TEST(DistanceCalculatorTests, EnforcesMinimumAndMaximumBondCutoffs) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("Si", {0.4, 0.0, 0.0}); // Distance = 0.4 Å (below min 0.8 Å)
  cell.addAtom("Si", {1.5, 0.0, 0.0}); // Distance = 1.5 Å (within [0.8, 2.0] Å)
  cell.addAtom("Si",
               {2.5, 0.0, 0.0}); // Distance = 2.5 Å (above max 2.0 Å, below global cutoff 3.0 Å)

  real_t const cutoff_sq = 9.0; // Global search cutoff = 3.0 Å

  // Bond cutoff window: [0.8 Å, 2.0 Å] -> [0.64 Å², 4.00 Å²]
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.64, .max_sq = 4.00}}};
  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(4);
  DistanceCalculationConfig const config{
      .r_max = 3.0,
      .r_bin_width = 0.02,
      .num_bins = 150,
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

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

TEST(DistanceCalculatorTests, LargeTriclinicCellSubcellPartitioning) {
  // Large non-orthogonal triclinic cell: 30x30x30 with alpha=75, beta=80, gamma=85 deg
  Cell cell({30.0, 30.0, 30.0, 75.0, 80.0, 85.0});

  // Atom 0 at origin
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // Atom 1 nearby in Cartesian space (distance = 2.0 Å)
  cell.addAtom("Si", {2.0, 0.0, 0.0});

  // Atom 2 placed far away in the center of the box (~15 Å away in a distant sub-cell)
  const auto &lattice = cell.latticeVectors();
  const auto center_pos = 0.5 * lattice[0] + 0.5 * lattice[1] + 0.5 * lattice[2];
  cell.addAtom("Si", center_pos);

  // Atom 3 placed near the opposite face to test sheared PBC boundary (wrapped distance ~ 1.5 Å)
  const auto near_pbc = lattice[0] - correlation::math::Vector3<real_t>{1.5, 0.0, 0.0};
  cell.addAtom("Si", near_pbc);

  real_t const cutoff = 3.5;
  real_t const cutoff_sq = cutoff * cutoff;
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.64, .max_sq = cutoff_sq}}};

  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(4);
  DistanceCalculationConfig const config{
      .r_max = cutoff,
      .r_bin_width = 0.02,
      .num_bins = static_cast<size_t>(std::ceil(cutoff / 0.02)),
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

  // Atom 0 <-> Atom 1 (2.0 Å) should be connected
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));

  // Atom 0 <-> Atom 2 (~15 Å, distant sub-cell) must NOT be connected
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{2}));
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{2}, AtomIndex{0}));

  // Atom 0 <-> Atom 3 (1.5 Å across PBC) should be connected
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{3}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{3}, AtomIndex{0}));
}

TEST(DistanceCalculatorTests, SmallCrystalCellMultiImageExpansion) {
  // 1-atom cubic crystal unit cell: a = 3.6 Å
  Cell cell({3.6, 0.0, 0.0}, {0.0, 3.6, 0.0}, {0.0, 0.0, 3.6});
  cell.addAtom("Cu", {0.0, 0.0, 0.0});

  real_t const r_max = 10.0;
  real_t const cutoff_sq = r_max * r_max;
  real_t const bin_width = 0.02;
  const auto num_bins = static_cast<size_t>(std::ceil(r_max / bin_width));

  // No bonds needed, only distance histogramming
  BondCutoffMatrix const empty_bonds = {{{.min_sq = 0.0, .max_sq = 0.0}}};
  RawHistogramTensor out_histograms;
  NeighborGraph out_graph(1);
  DistanceCalculationConfig const config{
      .r_max = r_max,
      .r_bin_width = bin_width,
      .num_bins = num_bins,
  };

  DistanceCalculator::compute(cell, cutoff_sq, empty_bonds, false, out_graph, &out_histograms,
                              config);

  ASSERT_EQ(out_histograms.size(), 1);
  ASSERT_EQ(out_histograms[0].size(), 1);

  // Shell 1: 6 neighbors at 3.6 Å (bin index 3.6 / 0.02 = 180, raw upper-triangle count = 3.0)
  auto const bin_shell1 = static_cast<size_t>(3.6 / bin_width);
  EXPECT_EQ(out_histograms[0][0][bin_shell1], 3.0);

  // Shell 2: 12 neighbors at 3.6 * sqrt(2) ≈ 5.09117 Å (raw upper-triangle count = 6.0)
  auto const bin_shell2 = static_cast<size_t>((3.6 * std::numbers::sqrt2) / bin_width);
  EXPECT_EQ(out_histograms[0][0][bin_shell2], 6.0);

  // Shell 3: 8 neighbors at 3.6 * sqrt(3) ≈ 6.23538 Å (raw upper-triangle count = 4.0)
  auto const bin_shell3 = static_cast<size_t>((3.6 * std::numbers::sqrt3) / bin_width);
  EXPECT_EQ(out_histograms[0][0][bin_shell3], 4.0);

  // Shell 4: 6 neighbors at 7.2 Å (raw upper-triangle count = 3.0)
  auto const bin_shell4 = static_cast<size_t>(7.2 / bin_width);
  EXPECT_EQ(out_histograms[0][0][bin_shell4], 3.0);
}

TEST(DistanceCalculatorTests, SmallCellFCCMultiImageExpansionMatchesSupercellUpTo20A) {
  real_t const lat_a = 3.6;
  real_t const r_max = 20.0;
  real_t const cutoff_sq = r_max * r_max;
  real_t const bin_width = 0.0217;
  const auto num_bins = static_cast<size_t>(std::ceil(r_max / bin_width));
  DistanceCalculationConfig const config{
      .r_max = r_max,
      .r_bin_width = bin_width,
      .num_bins = num_bins,
  };
  BondCutoffMatrix const empty_bonds = {{{.min_sq = 0.0, .max_sq = 0.0}}};

  // 1. Small unit cell: 1x1x1 FCC (4 atoms, cell dimension 3.6 Å)
  auto unit_cell = crystals::createFCCCell(lat_a, "Cu", 1, 1, 1);
  RawHistogramTensor unit_hists;
  NeighborGraph unit_graph(unit_cell.atoms().size());
  DistanceCalculator::compute(unit_cell, cutoff_sq, empty_bonds, false, unit_graph, &unit_hists,
                              config);

  // 2. Large supercell: 7x7x7 FCC (1372 atoms, cell dimension 25.2 Å > 20 Å)
  auto super_cell_7 = crystals::createFCCCell(lat_a, "Cu", 7, 7, 7);
  RawHistogramTensor super_hists_7;
  NeighborGraph super_graph_7(super_cell_7.atoms().size());
  DistanceCalculator::compute(super_cell_7, cutoff_sq, empty_bonds, false, super_graph_7,
                              &super_hists_7, config);

  // 3. Medium supercell: 5x5x5 FCC (500 atoms, cell dimension 18.0 Å)
  auto super_cell_5 = crystals::createFCCCell(lat_a, "Cu", 5, 5, 5);
  RawHistogramTensor super_hists_5;
  NeighborGraph super_graph_5(super_cell_5.atoms().size());
  DistanceCalculator::compute(super_cell_5, cutoff_sq, empty_bonds, false, super_graph_5,
                              &super_hists_5, config);

  const auto n_unit = static_cast<real_t>(unit_cell.atoms().size());
  const auto n_super_7 = static_cast<real_t>(super_cell_7.atoms().size());
  const auto n_super_5 = static_cast<real_t>(super_cell_5.atoms().size());

  ASSERT_EQ(unit_hists[0][0].size(), num_bins);
  ASSERT_EQ(super_hists_7[0][0].size(), num_bins);
  ASSERT_EQ(super_hists_5[0][0].size(), num_bins);

  real_t cum_unit = 0.0;
  real_t cum_super_7 = 0.0;
  real_t cum_super_5 = 0.0;
  for (size_t bin = 0; bin < num_bins; ++bin) {
    cum_unit += unit_hists[0][0][bin] / n_unit;
    cum_super_7 += super_hists_7[0][0][bin] / n_super_7;
    cum_super_5 += super_hists_5[0][0][bin] / n_super_5;
    EXPECT_NEAR(cum_unit, cum_super_7, 1e-3)
        << "Mismatch vs 7x7x7 at cumulative bin " << bin
        << " (r = " << (static_cast<real_t>(bin) + 1.0) * bin_width << " Å)";
    EXPECT_NEAR(cum_unit, cum_super_5, 1e-3)
        << "Mismatch vs 5x5x5 at cumulative bin " << bin
        << " (r = " << (static_cast<real_t>(bin) + 1.0) * bin_width << " Å)";
  }
}

TEST(DistanceCalculatorTests, SmallCellBCCMultiImageExpansionMatchesSupercellUpTo20A) {
  real_t const lat_a = 2.88;
  real_t const r_max = 20.0;
  real_t const cutoff_sq = r_max * r_max;
  real_t const bin_width = 0.03125;
  const auto num_bins = static_cast<size_t>(std::ceil(r_max / bin_width));
  DistanceCalculationConfig const config{
      .r_max = r_max,
      .r_bin_width = bin_width,
      .num_bins = num_bins,
  };
  BondCutoffMatrix const empty_bonds = {{{.min_sq = 0.0, .max_sq = 0.0}}};

  // 1. Small unit cell: 1x1x1 BCC (2 atoms, cell dimension 2.88 Å)
  auto unit_cell = crystals::createBCCCell(lat_a, "Fe", 1, 1, 1);
  RawHistogramTensor unit_hists;
  NeighborGraph unit_graph(unit_cell.atoms().size());
  DistanceCalculator::compute(unit_cell, cutoff_sq, empty_bonds, false, unit_graph, &unit_hists,
                              config);

  // 2. Large supercell: 7x7x7 BCC (686 atoms, cell dimension 20.16 Å > 20 Å)
  auto super_cell_7 = crystals::createBCCCell(lat_a, "Fe", 7, 7, 7);
  RawHistogramTensor super_hists_7;
  NeighborGraph super_graph_7(super_cell_7.atoms().size());
  DistanceCalculator::compute(super_cell_7, cutoff_sq, empty_bonds, false, super_graph_7,
                              &super_hists_7, config);

  // 3. Medium supercell: 5x5x5 BCC (250 atoms, cell dimension 14.4 Å)
  auto super_cell_5 = crystals::createBCCCell(lat_a, "Fe", 5, 5, 5);
  RawHistogramTensor super_hists_5;
  NeighborGraph super_graph_5(super_cell_5.atoms().size());
  DistanceCalculator::compute(super_cell_5, cutoff_sq, empty_bonds, false, super_graph_5,
                              &super_hists_5, config);

  const auto n_unit = static_cast<real_t>(unit_cell.atoms().size());
  const auto n_super_7 = static_cast<real_t>(super_cell_7.atoms().size());
  const auto n_super_5 = static_cast<real_t>(super_cell_5.atoms().size());

  ASSERT_EQ(unit_hists[0][0].size(), num_bins);
  ASSERT_EQ(super_hists_7[0][0].size(), num_bins);
  ASSERT_EQ(super_hists_5[0][0].size(), num_bins);

  real_t cum_unit = 0.0;
  real_t cum_super_7 = 0.0;
  real_t cum_super_5 = 0.0;
  for (size_t bin = 0; bin < num_bins; ++bin) {
    cum_unit += unit_hists[0][0][bin] / n_unit;
    cum_super_7 += super_hists_7[0][0][bin] / n_super_7;
    cum_super_5 += super_hists_5[0][0][bin] / n_super_5;
    EXPECT_NEAR(cum_unit, cum_super_7, 1e-3)
        << "Mismatch vs 7x7x7 at cumulative bin " << bin
        << " (r = " << (static_cast<real_t>(bin) + 1.0) * bin_width << " Å)";
    EXPECT_NEAR(cum_unit, cum_super_5, 1e-3)
        << "Mismatch vs 5x5x5 at cumulative bin " << bin
        << " (r = " << (static_cast<real_t>(bin) + 1.0) * bin_width << " Å)";
  }
}

TEST(DistanceCalculatorTests, Large40ACellSubcellPartitioningAndPBCWrapping) {
  // 40 Å x 40 Å x 40 Å cubic box
  Cell cell({40.0, 40.0, 40.0, 90.0, 90.0, 90.0});

  // Atom 0 at (0.5, 0.5, 0.5) in bin (0, 0, 0)
  cell.addAtom("Si", {0.5, 0.5, 0.5});

  // Atom 1 in same subcell at (1.5, 0.5, 0.5) (dist = 1.0 Å)
  cell.addAtom("Si", {1.5, 0.5, 0.5});

  // Atom 2 across X PBC boundary at (39.2, 0.5, 0.5) in bin (7, 0, 0) (wrapped dist = 1.3 Å)
  cell.addAtom("Si", {39.2, 0.5, 0.5});

  // Atom 3 across Y PBC boundary at (0.5, 39.0, 0.5) in bin (0, 7, 0) (wrapped dist = 1.5 Å)
  cell.addAtom("Si", {0.5, 39.0, 0.5});

  // Atom 4 across Z PBC boundary at (0.5, 0.5, 38.8) in bin (0, 0, 7) (wrapped dist = 1.7 Å)
  cell.addAtom("Si", {0.5, 0.5, 38.8});

  // Atom 5 across diagonal 3D corner PBC at (39.5, 39.5, 39.5) in bin (7, 7, 7) (wrapped dist =
  // sqrt(3) ≈ 1.732 Å)
  cell.addAtom("Si", {39.5, 39.5, 39.5});

  // Atom 6 deep in interior at (20.0, 20.0, 20.0) in bin (4, 4, 4) (~28 Å away from Atom 0)
  cell.addAtom("Si", {20.0, 20.0, 20.0});

  // Atom 7 near Atom 6 at (21.0, 20.0, 20.0) in bin (4, 4, 4) (dist = 1.0 Å)
  cell.addAtom("Si", {21.0, 20.0, 20.0});

  real_t const cutoff = 5.0;
  real_t const cutoff_sq = cutoff * cutoff;
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.25, .max_sq = cutoff_sq}}};

  NeighborGraph out_graph(8);
  RawHistogramTensor out_histograms;
  DistanceCalculationConfig const config{
      .r_max = cutoff,
      .r_bin_width = 0.05,
      .num_bins = static_cast<size_t>(std::ceil(cutoff / 0.05)),
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

  // Atom 0 <-> Atom 1 (1.0 Å within subcell)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));

  // Atom 0 <-> Atom 2 (1.3 Å across X PBC)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{2}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{2}, AtomIndex{0}));

  // Atom 0 <-> Atom 3 (1.5 Å across Y PBC)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{3}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{3}, AtomIndex{0}));

  // Atom 0 <-> Atom 4 (1.7 Å across Z PBC)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{4}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{4}, AtomIndex{0}));

  // Atom 0 <-> Atom 5 (1.732 Å across 3D corner PBC)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{5}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{5}, AtomIndex{0}));

  // Atom 0 must NOT connect to distant center atoms 6 or 7 (~28 Å away)
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{6}));
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{7}));

  // Atom 6 <-> Atom 7 (1.0 Å in interior subcell)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{6}, AtomIndex{7}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{7}, AtomIndex{6}));

  // Interior atoms must NOT connect to boundary atoms
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{6}, AtomIndex{2}));
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{6}, AtomIndex{5}));
}

TEST(DistanceCalculatorTests, Large40ATriclinicSubcellPartitioningAndPBCWrapping) {
  // 40 Å non-orthogonal triclinic box
  Cell cell({40.0, 40.0, 40.0, 75.0, 80.0, 85.0});
  const auto &lattice = cell.latticeVectors();

  // Atom 0 near origin
  cell.addAtom("Si", {0.2, 0.2, 0.2});

  // Atom 1 placed ~25 Å away in box center (distant subcell)
  const auto center_pos = 0.5 * lattice[0] + 0.5 * lattice[1] + 0.5 * lattice[2];
  cell.addAtom("Si", center_pos);

  // Atom 2 near A boundary face (wrapped distance ~1.3 Å)
  cell.addAtom("Si", lattice[0] - correlation::math::Vector3<real_t>{1.1, -0.2, -0.2});

  // Atom 3 near B boundary face (wrapped distance ~1.4 Å)
  cell.addAtom("Si", lattice[1] - correlation::math::Vector3<real_t>{-0.2, 1.2, -0.2});

  // Atom 4 near C boundary face (wrapped distance ~1.5 Å)
  cell.addAtom("Si", lattice[2] - correlation::math::Vector3<real_t>{-0.2, -0.2, 1.3});

  real_t const cutoff = 4.0;
  real_t const cutoff_sq = cutoff * cutoff;
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.25, .max_sq = cutoff_sq}}};

  NeighborGraph out_graph(5);
  RawHistogramTensor out_histograms;
  DistanceCalculationConfig const config{
      .r_max = cutoff,
      .r_bin_width = 0.05,
      .num_bins = static_cast<size_t>(std::ceil(cutoff / 0.05)),
  };

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms,
                              config);

  // Atom 0 <-> Atom 1 (distant center subcell) must NOT be connected
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{0}, AtomIndex{1}));
  EXPECT_FALSE(out_graph.areConnected(AtomIndex{1}, AtomIndex{0}));

  // Atom 0 <-> Atom 2 (across A PBC face)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{2}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{2}, AtomIndex{0}));

  // Atom 0 <-> Atom 3 (across B PBC face)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{3}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{3}, AtomIndex{0}));

  // Atom 0 <-> Atom 4 (across C PBC face)
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{0}, AtomIndex{4}));
  EXPECT_TRUE(out_graph.areConnected(AtomIndex{4}, AtomIndex{0}));
}

// --- Brute-Force Reference Equivalence Tests ---

namespace {

struct BruteForceReferenceData {
  std::vector<real_t> distances;
  RawHistogramTensor histograms;
  size_t total_directed_bonds = 0;
};

std::vector<math::Vector3<real_t>> getWrappedPositions(const Cell &cell) {
  const auto &atoms = cell.atoms();
  const auto &lattice = cell.latticeVectors();
  const auto inv = cell.inverseLatticeVectors();
  std::vector<math::Vector3<real_t>> wrapped;
  wrapped.reserve(atoms.size());

  for (const auto &atom : atoms) {
    const math::Vector3<real_t> frac = inv * atom.position();
    const real_t frac_x = frac.x() - std::floor(frac.x());
    const real_t frac_y = frac.y() - std::floor(frac.y());
    const real_t frac_z = frac.z() - std::floor(frac.z());
    wrapped.push_back(frac_x * lattice[0] + frac_y * lattice[1] + frac_z * lattice[2]);
  }
  return wrapped;
}

std::tuple<int, int, int> computeMaxShifts(const Cell &cell, real_t cutoff) {
  const real_t vol = std::abs(cell.volume());
  const auto &lat = cell.latticeVectors();
  const real_t width_x = vol / math::norm(math::cross(lat[1], lat[2]));
  const real_t width_y = vol / math::norm(math::cross(lat[2], lat[0]));
  const real_t width_z = vol / math::norm(math::cross(lat[0], lat[1]));

  const int shift_x = std::max(1, static_cast<int>(std::ceil(cutoff / width_x)) + 1);
  const int shift_y = std::max(1, static_cast<int>(std::ceil(cutoff / width_y)) + 1);
  const int shift_z = std::max(1, static_cast<int>(std::ceil(cutoff / width_z)) + 1);
  return {shift_x, shift_y, shift_z};
}

void processCandidatePair(size_t first_idx, size_t second_idx, const math::Vector3<real_t> &r_ij,
                          real_t cutoff_sq, const Cell &cell, const BondCutoffMatrix &bond_cutoffs,
                          const DistanceCalculationConfig &config,
                          BruteForceReferenceData &ref_data) {
  const real_t d_sq = math::norm_sq(r_ij);
  if (d_sq >= cutoff_sq) {
    return;
  }

  const real_t dist = std::sqrt(d_sq);
  ref_data.distances.push_back(dist);

  const int type_a = cell.atoms()[first_idx].elementId();
  const int type_b = cell.atoms()[second_idx].elementId();

  if (config.num_bins > 0 && config.r_bin_width > 0.0 && dist < config.r_max) {
    const auto bin = static_cast<size_t>(dist / config.r_bin_width);
    if (bin < config.num_bins) {
      ref_data.histograms[type_a][type_b][bin] += 1.0;
      if (first_idx != second_idx && type_a != type_b) {
        ref_data.histograms[type_b][type_a][bin] += 1.0;
      }
    }
  }

  if (static_cast<size_t>(type_a) < bond_cutoffs.size() &&
      static_cast<size_t>(type_b) < bond_cutoffs[type_a].size()) {
    const auto &range = bond_cutoffs[type_a][type_b];
    if (d_sq >= range.min_sq && d_sq <= range.max_sq) {
      ref_data.total_directed_bonds += 2;
    }
  }
}

void evaluateShiftedPairs(size_t atom_idx_a, size_t atom_idx_b, const math::Vector3<real_t> &pos_a,
                          const math::Vector3<real_t> &pos_b, const math::Vector3<real_t> &shift,
                          real_t cutoff_sq, const Cell &cell, const BondCutoffMatrix &bond_cutoffs,
                          const DistanceCalculationConfig &config,
                          BruteForceReferenceData &ref_data) {
  if (atom_idx_a == atom_idx_b) {
    processCandidatePair(atom_idx_a, atom_idx_b, shift, cutoff_sq, cell, bond_cutoffs, config,
                         ref_data);
  } else {
    processCandidatePair(atom_idx_a, atom_idx_b, (pos_b + shift) - pos_a, cutoff_sq, cell,
                         bond_cutoffs, config, ref_data);
    processCandidatePair(atom_idx_b, atom_idx_a, (pos_a + shift) - pos_b, cutoff_sq, cell,
                         bond_cutoffs, config, ref_data);
  }
}

void evaluateSingleShift(int shift_dx, int shift_dy, int shift_dz, const Cell &cell,
                         const std::vector<math::Vector3<real_t>> &wrapped, real_t cutoff_sq,
                         const BondCutoffMatrix &bond_cutoffs,
                         bool ignore_periodic_self_interactions,
                         const DistanceCalculationConfig &config,
                         BruteForceReferenceData &ref_data) {
  const bool is_positive = (shift_dx > 0) || (shift_dx == 0 && shift_dy > 0) ||
                           (shift_dx == 0 && shift_dy == 0 && shift_dz > 0);
  if (!is_positive) {
    return;
  }

  const auto &lat = cell.latticeVectors();
  const math::Vector3<real_t> shift = static_cast<real_t>(shift_dx) * lat[0] +
                                      static_cast<real_t>(shift_dy) * lat[1] +
                                      static_cast<real_t>(shift_dz) * lat[2];

  const size_t n_atoms = cell.atomCount();
  for (size_t atom_idx_a = 0; atom_idx_a < n_atoms; ++atom_idx_a) {
    const size_t start_b = ignore_periodic_self_interactions ? (atom_idx_a + 1) : atom_idx_a;
    for (size_t atom_idx_b = start_b; atom_idx_b < n_atoms; ++atom_idx_b) {
      evaluateShiftedPairs(atom_idx_a, atom_idx_b, wrapped[atom_idx_a], wrapped[atom_idx_b], shift,
                           cutoff_sq, cell, bond_cutoffs, config, ref_data);
    }
  }
}

void evaluateAllPeriodicShifts(const Cell &cell, int max_dx, int max_dy, int max_dz,
                               const std::vector<math::Vector3<real_t>> &wrapped, real_t cutoff_sq,
                               const BondCutoffMatrix &bond_cutoffs,
                               bool ignore_periodic_self_interactions,
                               const DistanceCalculationConfig &config,
                               BruteForceReferenceData &ref_data) {
  for (int shift_dx = -max_dx; shift_dx <= max_dx; ++shift_dx) {
    for (int shift_dy = -max_dy; shift_dy <= max_dy; ++shift_dy) {
      for (int shift_dz = -max_dz; shift_dz <= max_dz; ++shift_dz) {
        evaluateSingleShift(shift_dx, shift_dy, shift_dz, cell, wrapped, cutoff_sq, bond_cutoffs,
                            ignore_periodic_self_interactions, config, ref_data);
      }
    }
  }
}

BruteForceReferenceData computeBruteForceReference(const Cell &cell, real_t cutoff_sq,
                                                   const BondCutoffMatrix &bond_cutoffs,
                                                   bool ignore_periodic_self_interactions,
                                                   DistanceCalculationConfig config) {
  BruteForceReferenceData ref_data;
  const size_t num_elements = cell.elements().size();
  if (config.num_bins > 0) {
    ref_data.histograms.assign(
        num_elements,
        std::vector<std::vector<real_t>>(num_elements, std::vector<real_t>(config.num_bins, 0.0)));
  }

  const auto wrapped = getWrappedPositions(cell);
  const auto [max_dx, max_dy, max_dz] = computeMaxShifts(cell, std::sqrt(cutoff_sq));
  const size_t n_atoms = cell.atomCount();

  // Primary cell pairs (dx = dy = dz = 0)
  for (size_t atom_idx_a = 0; atom_idx_a < n_atoms; ++atom_idx_a) {
    for (size_t atom_idx_b = atom_idx_a + 1; atom_idx_b < n_atoms; ++atom_idx_b) {
      processCandidatePair(atom_idx_a, atom_idx_b, wrapped[atom_idx_b] - wrapped[atom_idx_a],
                           cutoff_sq, cell, bond_cutoffs, config, ref_data);
    }
  }

  evaluateAllPeriodicShifts(cell, max_dx, max_dy, max_dz, wrapped, cutoff_sq, bond_cutoffs,
                            ignore_periodic_self_interactions, config, ref_data);

  std::ranges::sort(ref_data.distances);
  return ref_data;
}

void verifyBonds(const Cell &cell, const NeighborGraph &actual_graph,
                 const BruteForceReferenceData &ref_data, std::vector<real_t> &actual_distances) {
  size_t actual_bonds = 0;
  for (size_t atom_idx = 0; atom_idx < cell.atomCount(); ++atom_idx) {
    const auto &nbrs = actual_graph.getNeighbors(atom_idx);
    actual_bonds += nbrs.size();
    for (const auto &nbr : nbrs) {
      if (atom_idx < nbr.index) {
        actual_distances.push_back(nbr.distance);
      }
    }
  }
  EXPECT_EQ(actual_bonds, ref_data.total_directed_bonds);
}

void verifyHistogramBins(const std::vector<real_t> &actual_bins,
                         const std::vector<real_t> &ref_bins, size_t elem_a, size_t elem_b) {
  ASSERT_EQ(actual_bins.size(), ref_bins.size());
  for (size_t bin = 0; bin < actual_bins.size(); ++bin) {
    EXPECT_EQ(actual_bins[bin], ref_bins[bin])
        << "Mismatch at [" << elem_a << "][" << elem_b << "][" << bin << "]";
  }
}

void verifyHistograms(const RawHistogramTensor &actual_histograms,
                      const RawHistogramTensor &ref_histograms) {
  ASSERT_EQ(actual_histograms.size(), ref_histograms.size());
  for (size_t elem_a = 0; elem_a < actual_histograms.size(); ++elem_a) {
    ASSERT_EQ(actual_histograms[elem_a].size(), ref_histograms[elem_a].size());
    for (size_t elem_b = 0; elem_b < actual_histograms[elem_a].size(); ++elem_b) {
      verifyHistogramBins(actual_histograms[elem_a][elem_b], ref_histograms[elem_a][elem_b], elem_a,
                          elem_b);
    }
  }
}

void verifyDistances(std::vector<real_t> &actual_distances,
                     const std::vector<real_t> &ref_distances) {
  std::ranges::sort(actual_distances);
  ASSERT_EQ(actual_distances.size(), ref_distances.size());
  for (size_t idx = 0; idx < actual_distances.size(); ++idx) {
    const real_t expected = ref_distances[idx];
    const real_t actual = actual_distances[idx];
    const real_t rel_eps =
        correlation::is_single_precision ? static_cast<real_t>(1e-5) : static_cast<real_t>(1e-12);
    const real_t tol = std::max(rel_eps, rel_eps * expected);
    EXPECT_NEAR(actual, expected, tol) << "Distance mismatch at index " << idx;
  }
}

void verifyDistanceEquivalence(const Cell &cell, real_t cutoff,
                               const BondCutoffMatrix &bond_cutoffs,
                               const DistanceCalculationConfig &config) {
  const real_t cutoff_sq = cutoff * cutoff;
  BondCutoffMatrix effective_cutoffs = bond_cutoffs;
  for (auto &row : effective_cutoffs) {
    for (auto &range : row) {
      range.max_sq = std::min(range.max_sq, cutoff_sq);
    }
  }
  NeighborGraph actual_graph(cell.atomCount());
  RawHistogramTensor actual_histograms;

  DistanceCalculator::compute(cell, cutoff_sq, effective_cutoffs, true, actual_graph,
                              &actual_histograms, config);

  const auto ref = computeBruteForceReference(cell, cutoff_sq, effective_cutoffs, true, config);

  std::vector<real_t> actual_distances;
  verifyBonds(cell, actual_graph, ref, actual_distances);
  verifyHistograms(actual_histograms, ref.histograms);
  verifyDistances(actual_distances, ref.distances);
}

struct RandomAmorphousCellConfig {
  size_t num_atoms = 0;
  real_t box_len = 0.0;
  uint32_t seed = 0;
};

Cell makeRandomAmorphousCell(const RandomAmorphousCellConfig &random_cfg) {
  Cell cell({random_cfg.box_len, 0.0, 0.0}, {0.0, random_cfg.box_len, 0.0},
            {0.0, 0.0, random_cfg.box_len});
  std::mt19937 gen(random_cfg.seed);
  std::uniform_real_distribution<real_t> dist(static_cast<real_t>(0.1),
                                              random_cfg.box_len - static_cast<real_t>(0.1));
  for (size_t idx = 0; idx < random_cfg.num_atoms; ++idx) {
    const real_t coord_x = dist(gen);
    const real_t coord_y = dist(gen);
    const real_t coord_z = dist(gen);
    cell.addAtom("Si", math::Vector3<real_t>{coord_x, coord_y, coord_z});
  }
  return cell;
}

} // namespace

TEST(DistanceCalculatorTests, BruteForceEquivalenceTinySingleAtomCell) {
  Cell cell({3.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {0.0, 0.0, 3.0});
  cell.addAtom("Ar", {1.5, 1.5, 1.5});
  const BondCutoffMatrix bonds = {{{.min_sq = 0.0, .max_sq = 100.0}}};
  const DistanceCalculationConfig config{.r_max = 8.0, .r_bin_width = 0.1, .num_bins = 80};

  for (const real_t cutoff :
       {static_cast<real_t>(2.5), static_cast<real_t>(5.0), static_cast<real_t>(7.5)}) {
    verifyDistanceEquivalence(cell, cutoff, bonds, config);
  }
}

TEST(DistanceCalculatorTests, BruteForceEquivalenceSkewedTriclinicCell) {
  Cell cell({8.0, 9.0, 10.0, 65.0, 75.0, 85.0});
  cell.addAtom("Si", {0.5, 0.5, 0.5});
  cell.addAtom("Si", {2.0, 3.0, 1.5});
  cell.addAtom("O", {4.0, 2.0, 5.0});
  cell.addAtom("O", {6.0, 7.0, 8.0});

  const BondCutoffMatrix bonds = {
      {{.min_sq = 0.0, .max_sq = 100.0}, {.min_sq = 0.0, .max_sq = 100.0}},
      {{.min_sq = 0.0, .max_sq = 100.0}, {.min_sq = 0.0, .max_sq = 100.0}},
  };
  const DistanceCalculationConfig config{.r_max = 6.0, .r_bin_width = 0.1, .num_bins = 60};

  for (const real_t cutoff : {static_cast<real_t>(3.0), static_cast<real_t>(5.5)}) {
    verifyDistanceEquivalence(cell, cutoff, bonds, config);
  }
}

TEST(DistanceCalculatorTests, BruteForceEquivalenceThinSlabCell) {
  // One axis is thin (3.0 Å), other two axes are 16.0 Å
  Cell cell({16.0, 0.0, 0.0}, {0.0, 16.0, 0.0}, {0.0, 0.0, 3.0});
  cell.addAtom("Cu", {2.0, 2.0, 1.5});
  cell.addAtom("Cu", {5.0, 8.0, 1.0});
  cell.addAtom("Cu", {14.0, 14.0, 2.0});

  const BondCutoffMatrix bonds = {{{.min_sq = 0.0, .max_sq = 100.0}}};
  const DistanceCalculationConfig config{.r_max = 7.0, .r_bin_width = 0.1, .num_bins = 70};

  for (const real_t cutoff : {static_cast<real_t>(4.0), static_cast<real_t>(6.5)}) {
    verifyDistanceEquivalence(cell, cutoff, bonds, config);
  }
}

TEST(DistanceCalculatorTests, BruteForceEquivalenceBoundaryWidthCell) {
  // Cubic cell width = 6.0 Å. Test at cutoff = 6.0 (width == cutoff) and cutoff = 3.0 (width ==
  // 2*cutoff)
  Cell cell({6.0, 0.0, 0.0}, {0.0, 6.0, 0.0}, {0.0, 0.0, 6.0});
  cell.addAtom("Fe", {0.5, 0.5, 0.5});
  cell.addAtom("Fe", {3.0, 3.0, 3.0});
  cell.addAtom("Fe", {5.5, 2.0, 4.0});

  const BondCutoffMatrix bonds = {{{.min_sq = 0.0, .max_sq = 100.0}}};
  const DistanceCalculationConfig config{.r_max = 6.0, .r_bin_width = 0.05, .num_bins = 120};

  for (const real_t cutoff : {static_cast<real_t>(3.0), static_cast<real_t>(6.0)}) {
    verifyDistanceEquivalence(cell, cutoff, bonds, config);
  }
}

TEST(DistanceCalculatorTests, BruteForceEquivalenceDiamondAndFCCSupercells) {
  const auto diamond = crystals::createDiamondCell(5.43, "Si", 2, 2, 2);
  const auto fcc = crystals::createFCCCell(3.61, "Cu", 2, 2, 2);

  const BondCutoffMatrix si_bonds = {{{.min_sq = 0.0, .max_sq = 25.0}}};
  const DistanceCalculationConfig si_config{.r_max = 4.5, .r_bin_width = 0.05, .num_bins = 90};
  verifyDistanceEquivalence(diamond, static_cast<real_t>(4.0), si_bonds, si_config);

  const BondCutoffMatrix cu_bonds = {{{.min_sq = 0.0, .max_sq = 25.0}}};
  const DistanceCalculationConfig cu_config{.r_max = 4.5, .r_bin_width = 0.05, .num_bins = 90};
  verifyDistanceEquivalence(fcc, static_cast<real_t>(3.5), cu_bonds, cu_config);
}

TEST(DistanceCalculatorTests, BruteForceEquivalenceRandomAmorphousCell) {
  const auto cell = makeRandomAmorphousCell(RandomAmorphousCellConfig{
      .num_atoms = 20,
      .box_len = static_cast<real_t>(10.0),
      .seed = 42,
  });
  const BondCutoffMatrix bonds = {{{.min_sq = 0.0, .max_sq = 50.0}}};
  const DistanceCalculationConfig config{.r_max = 5.0, .r_bin_width = 0.05, .num_bins = 100};
  for (const real_t cutoff : {static_cast<real_t>(2.5), static_cast<real_t>(4.8)}) {
    verifyDistanceEquivalence(cell, cutoff, bonds, config);
  }
}

TEST(DistanceCalculatorTests, DeterministicNeighborGraphAcrossThreadCounts) {
  const auto cell = makeRandomAmorphousCell(RandomAmorphousCellConfig{
      .num_atoms = 150,
      .box_len = static_cast<real_t>(18.0),
      .seed = 12345,
  });
  const real_t cutoff = 4.5;
  const real_t cutoff_sq = cutoff * cutoff;
  const BondCutoffMatrix bonds = {{{.min_sq = 0.0, .max_sq = cutoff_sq}}};
  const DistanceCalculationConfig config{.r_max = cutoff, .r_bin_width = 0.05, .num_bins = 90};

  NeighborGraph graph_1thread(cell.atomCount());
  {
    tbb::global_control global_control_1thread(tbb::global_control::max_allowed_parallelism, 1);
    DistanceCalculator::compute(cell, cutoff_sq, bonds, true, graph_1thread, nullptr, config);
  }

  NeighborGraph graph_multithread(cell.atomCount());
  {
    tbb::global_control global_control_multithread(tbb::global_control::max_allowed_parallelism, 8);
    DistanceCalculator::compute(cell, cutoff_sq, bonds, true, graph_multithread, nullptr, config);
  }

  ASSERT_EQ(graph_1thread.nodeCount(), graph_multithread.nodeCount());
  ASSERT_EQ(graph_1thread.edgeCount(), graph_multithread.edgeCount());

  for (size_t i = 0; i < cell.atomCount(); ++i) {
    auto nbrs1 = graph_1thread.getNeighbors(i);
    auto nbrs2 = graph_multithread.getNeighbors(i);
    ASSERT_EQ(nbrs1.size(), nbrs2.size()) << "Mismatch in neighbor count for atom " << i;

    for (size_t k = 0; k < nbrs1.size(); ++k) {
      EXPECT_EQ(nbrs1[k].index, nbrs2[k].index)
          << "Index mismatch for atom " << i << " at neighbor " << k;
      EXPECT_DOUBLE_EQ(nbrs1[k].distance, nbrs2[k].distance)
          << "Distance mismatch for atom " << i << " at neighbor " << k;
      EXPECT_DOUBLE_EQ(nbrs1[k].r_ij.x(), nbrs2[k].r_ij.x());
      EXPECT_DOUBLE_EQ(nbrs1[k].r_ij.y(), nbrs2[k].r_ij.y());
      EXPECT_DOUBLE_EQ(nbrs1[k].r_ij.z(), nbrs2[k].r_ij.z());
    }

    auto indices1 = graph_1thread.getNeighborIndices(i);
    auto indices2 = graph_multithread.getNeighborIndices(i);
    ASSERT_EQ(indices1.size(), indices2.size());
    for (size_t k = 0; k < indices1.size(); ++k) {
      EXPECT_EQ(indices1[k], indices2[k]);
    }
  }
}

TEST(DistanceCalculatorTests, ForEachPairWithinVisitsExpectedPairs) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});
  cell.addAtom("O", {2.0, 0.0, 0.0});
  cell.addAtom("O", {0.0, 5.0, 0.0});

  std::atomic<size_t> pair_count{0};
  DistanceCalculator::forEachPairWithin(cell, 3.0, true,
                                        [&](size_t atom_i, size_t atom_j, real_t d_sq,
                                            const correlation::math::Vector3<real_t> &r_ij) {
                                          pair_count.fetch_add(1, std::memory_order_relaxed);
                                          EXPECT_NE(atom_i, atom_j);
                                          EXPECT_LE(d_sq, 9.0);
                                          EXPECT_NEAR(std::sqrt(d_sq), 2.0, 1e-5);
                                          EXPECT_NEAR(std::abs(r_ij.x()), 2.0, 1e-5);
                                        });

  EXPECT_EQ(pair_count.load(), 1);

  // With r_max = 6.0, (0, 1) at 2.0 and (0, 2) at 5.0 are within cutoff
  pair_count.store(0);
  DistanceCalculator::forEachPairWithin(cell, 6.0, true,
                                        [&](size_t atom_i, size_t atom_j, real_t d_sq,
                                            const correlation::math::Vector3<real_t> &r_ij) {
                                          pair_count.fetch_add(1, std::memory_order_relaxed);
                                          EXPECT_NE(atom_i, atom_j);
                                          EXPECT_LE(d_sq, 36.0);
                                          EXPECT_GT(d_sq, 0.0);
                                          EXPECT_TRUE(r_ij.x() != 0.0 || r_ij.y() != 0.0 ||
                                                      r_ij.z() != 0.0);
                                        });

  EXPECT_GE(pair_count.load(), 2);
}

TEST(DistanceCalculatorTests, ForEachPairWithinHandlesEmptyAndSingleAtomCell) {
  const Cell empty_cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  size_t count = 0;
  DistanceCalculator::forEachPairWithin(empty_cell, 5.0, true,
                                        [&](size_t, size_t, real_t, const auto &) { ++count; });
  EXPECT_EQ(count, 0);

  Cell single_atom_cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  single_atom_cell.addAtom("Si", {5.0, 5.0, 5.0});
  DistanceCalculator::forEachPairWithin(single_atom_cell, 5.0, true,
                                        [&](size_t, size_t, real_t, const auto &) { ++count; });
  EXPECT_EQ(count, 0);

  // Negative cutoff
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {1.0, 0.0, 0.0});
  cell.addAtom("Si", {2.0, 0.0, 0.0});
  DistanceCalculator::forEachPairWithin(cell, -1.0, true,
                                        [&](size_t, size_t, real_t, const auto &) { ++count; });
  EXPECT_EQ(count, 0);
}

} // namespace correlation::testing
