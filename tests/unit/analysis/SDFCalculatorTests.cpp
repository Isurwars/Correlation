#include <gtest/gtest.h>
#include "calculators/SDFCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

using namespace correlation::calculators;
using namespace correlation::analysis;

TEST(SDFCalculatorTests, CalculateSDF) {
    correlation::core::Cell cell;
    cell.setLatticeParameters({10.0, 10.0, 10.0, 90.0, 90.0, 90.0});
    cell.addAtom("Ar", {5.0, 5.0, 5.0});

    DistributionFunctions df(cell);
    AnalysisSettings settings;
    settings.r_bin_width = 1.0;

    SDFCalculator calc;
    calc.calculateFrame(df, settings);

    const auto& hist = df.getHistogram("SDF");
    EXPECT_EQ(hist.bins.size(), 1000); // 10x10x10 grid
    EXPECT_TRUE(hist.partials.count("Ar") > 0);
    EXPECT_TRUE(hist.partials.count("Total") > 0);
}
