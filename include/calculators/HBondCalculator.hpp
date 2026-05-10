/**
 * @file HBondCalculator.hpp
 * @brief Hydrogen Bond Analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#pragma once

#include "BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @class HBondCalculator
 * @brief Computes hydrogen bond statistics based on geometric criteria.
 */
class HBondCalculator : public BaseCalculator {
public:
    std::string getName() const override { return "Hydrogen Bond"; }
    std::string getShortName() const override { return "HBond"; }
    std::string getGroup() const override { return "Structural"; }
    std::string getDescription() const override {
        return "Computes hydrogen bond statistics (counts and distribution).";
    }

    bool isFrameCalculator() const override { return true; }
    bool isTrajectoryCalculator() const override { return false; }

    void calculateFrame(
        correlation::analysis::DistributionFunctions &df,
        const correlation::analysis::AnalysisSettings &settings) const override;

    /**
     * @brief Performs Hydrogen Bond analysis.
     * @param cell The simulation cell.
     * @param neighbors The structural analyzer.
     * @return A histogram of H-bond counts.
     */
    static correlation::analysis::Histogram calculate(
        const correlation::core::Cell& cell,
        const correlation::analysis::StructureAnalyzer* neighbors);
};

} // namespace correlation::calculators
