#include "analysis/AnalysisTypes.hpp"
#include "calculators/DistanceCalculator.hpp"
#include "core/Cell.hpp"
#include "core/NeighborGraph.hpp"

#include "../../CrystalTestHelper.hpp"

#include <gtest/gtest.h>
#include <numbers>
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
      {{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}}  // O-Si, O-O
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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 4.0}}};

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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 25.0}}};

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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 9.0}}};

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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}},
                                         {{.min_sq = 0.36, .max_sq = 4.0}, {.min_sq = 0.36, .max_sq = 4.0}}};

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

  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.36, .max_sq = 4.0}}};
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, -1.0, bond_cutoffs, true, out_graph), std::invalid_argument);
  EXPECT_THROW(DistanceCalculator::compute(cell, 0.0, bond_cutoffs, true, out_graph), std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenMinBondCutoffExceedsMax) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  // min_sq (4.0) > max_sq (2.25)
  BondCutoffMatrix const invalid_cutoffs = {{{.min_sq = 4.0, .max_sq = 2.25}}};
  NeighborGraph out_graph(1);

  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, invalid_cutoffs, true, out_graph), std::invalid_argument);
}

TEST(DistanceCalculatorTests, ThrowsWhenCutoffBoundsAreNegative) {
  Cell cell({10.0, 0.0, 0.0}, {0.0, 10.0, 0.0}, {0.0, 0.0, 10.0});
  cell.addAtom("Si", {0.0, 0.0, 0.0});

  NeighborGraph out_graph(1);

  BondCutoffMatrix const neg_min = {{{.min_sq = -0.5, .max_sq = 4.0}}};
  EXPECT_THROW(DistanceCalculator::compute(cell, 9.0, neg_min, true, out_graph), std::invalid_argument);

  BondCutoffMatrix const neg_max = {{{.min_sq = 0.5, .max_sq = -4.0}}};
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
  BondCutoffMatrix const bond_cutoffs = {{{.min_sq = 0.64, .max_sq = 4.00}}};
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

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

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

  DistanceCalculator::compute(cell, cutoff_sq, empty_bonds, false, out_graph, &out_histograms, config);

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
  DistanceCalculator::compute(unit_cell, cutoff_sq, empty_bonds, false, unit_graph, &unit_hists, config);

  // 2. Large supercell: 7x7x7 FCC (1372 atoms, cell dimension 25.2 Å > 20 Å)
  auto super_cell_7 = crystals::createFCCCell(lat_a, "Cu", 7, 7, 7);
  RawHistogramTensor super_hists_7;
  NeighborGraph super_graph_7(super_cell_7.atoms().size());
  DistanceCalculator::compute(super_cell_7, cutoff_sq, empty_bonds, false, super_graph_7, &super_hists_7, config);

  // 3. Medium supercell: 5x5x5 FCC (500 atoms, cell dimension 18.0 Å)
  auto super_cell_5 = crystals::createFCCCell(lat_a, "Cu", 5, 5, 5);
  RawHistogramTensor super_hists_5;
  NeighborGraph super_graph_5(super_cell_5.atoms().size());
  DistanceCalculator::compute(super_cell_5, cutoff_sq, empty_bonds, false, super_graph_5, &super_hists_5, config);

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
    EXPECT_NEAR(cum_unit, cum_super_7, 1e-3) << "Mismatch vs 7x7x7 at cumulative bin " << bin
                                             << " (r = " << (static_cast<real_t>(bin) + 1.0) * bin_width << " Å)";
    EXPECT_NEAR(cum_unit, cum_super_5, 1e-3) << "Mismatch vs 5x5x5 at cumulative bin " << bin
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
  DistanceCalculator::compute(unit_cell, cutoff_sq, empty_bonds, false, unit_graph, &unit_hists, config);

  // 2. Large supercell: 7x7x7 BCC (686 atoms, cell dimension 20.16 Å > 20 Å)
  auto super_cell_7 = crystals::createBCCCell(lat_a, "Fe", 7, 7, 7);
  RawHistogramTensor super_hists_7;
  NeighborGraph super_graph_7(super_cell_7.atoms().size());
  DistanceCalculator::compute(super_cell_7, cutoff_sq, empty_bonds, false, super_graph_7, &super_hists_7, config);

  // 3. Medium supercell: 5x5x5 BCC (250 atoms, cell dimension 14.4 Å)
  auto super_cell_5 = crystals::createBCCCell(lat_a, "Fe", 5, 5, 5);
  RawHistogramTensor super_hists_5;
  NeighborGraph super_graph_5(super_cell_5.atoms().size());
  DistanceCalculator::compute(super_cell_5, cutoff_sq, empty_bonds, false, super_graph_5, &super_hists_5, config);

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
    EXPECT_NEAR(cum_unit, cum_super_7, 1e-3) << "Mismatch vs 7x7x7 at cumulative bin " << bin
                                             << " (r = " << (static_cast<real_t>(bin) + 1.0) * bin_width << " Å)";
    EXPECT_NEAR(cum_unit, cum_super_5, 1e-3) << "Mismatch vs 5x5x5 at cumulative bin " << bin
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

  // Atom 5 across diagonal 3D corner PBC at (39.5, 39.5, 39.5) in bin (7, 7, 7) (wrapped dist = sqrt(3) ≈ 1.732 Å)
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

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

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

  DistanceCalculator::compute(cell, cutoff_sq, bond_cutoffs, true, out_graph, &out_histograms, config);

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

} // namespace correlation::testing
