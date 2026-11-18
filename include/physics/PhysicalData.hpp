/**
 * @file PhysicalData.hpp
 * @brief Database of physical properties for chemical elements (radii, masses,
 * form factors).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

/**
 * @namespace correlation::physics
 * @brief Provides structures and functions for handling atomic and material
 * physical properties.
 */
namespace correlation::physics {

namespace detail {
/**
 * @struct ElementData
 * @brief Internal structure to hold element-specific physical data.
 *
 * This data is used for calculating density, scattering intensities, and
 * neighbor interactions.
 */
struct ElementData {
  std::string_view symbol; ///< Chemical symbol (e.g., "Si").
  real_t radius;           ///< Covalent radius in Angstroms.
  real_t mass;             ///< Atomic mass in atomic mass units (Da).

  /**
   * @brief 9 coefficients (a1, b1, a2, b2, a3, b3, a4, b4, c) for the
   * Cromer-Mann equation used to calculate atomic form factors f(q).
   */
  std::array<real_t, 9> form_factors;
};

/**
 * @brief Searches for an element by its symbol in the static database.
 * @param symbol The chemical symbol to look for.
 * @return Pointer to ElementData if found, nullptr otherwise.
 */
constexpr const ElementData *find(std::string_view symbol) noexcept;
} // namespace detail

/**
 * @brief Retrieves the covalent radius of a given chemical element.
 *
 * @param symbol The atomic symbol of the element (e.g., "C", "O", "Fe").
 * @return The covalent radius in Angstroms.
 * @throws std::out_of_range if the element symbol is not found in the database.
 */
inline real_t getCovalentRadius(std::string_view symbol) {
  const auto *data = detail::find(symbol);
  if (data == nullptr) {
    throw std::out_of_range(std::string("Covalent radius not found for element: ") +
                            std::string(symbol));
  }
  return data->radius;
}

/**
 * @brief Retrieves the atomic mass of a given chemical element without throwing.
 *
 * @param[in] symbol The atomic symbol of the element.
 * @return The atomic mass in atomic mass units (Da or g/mol), or std::nullopt if unknown.
 */
[[nodiscard]] constexpr std::optional<real_t> findAtomicMass(std::string_view symbol) noexcept {
  const auto *data = detail::find(symbol);
  if (data == nullptr) {
    return std::nullopt;
  }
  return data->mass;
}

/**
 * @brief Retrieves the atomic mass of a given chemical element.
 *
 * @param symbol The atomic symbol of the element.
 * @return The atomic mass in atomic mass units (Da or g/mol).
 * @throws std::out_of_range if the element symbol is not found in the database.
 */
inline real_t getAtomicMass(std::string_view symbol) {
  const auto *data = detail::find(symbol);
  if (data == nullptr) {
    throw std::out_of_range(std::string("Atomic mass not found for element: ") +
                            std::string(symbol));
  }
  return data->mass;
}

/**
 * @brief Retrieves the atomic form factors for a given chemical element.
 *
 * The array contains 9 coefficients (a1, b1, a2, b2, a3, b3, a4, b4, c) based
 * on the Cromer-Mann equation.
 *
 * @param symbol The atomic symbol of the element.
 * @return A 9-element array containing the form factor coefficients.
 * @throws std::out_of_range if the element symbol is not found in the database.
 */
inline const std::array<real_t, 9> &getAtomicFormFactors(std::string_view symbol) {
  const auto *data = detail::find(symbol);
  if (data == nullptr) {
    throw std::out_of_range(std::string("Atomic form factor not found for element: ") +
                            std::string(symbol));
  }
  return data->form_factors;
}

namespace detail {
inline constexpr std::array<ElementData, 98> DATABASE = {
    {
        {
            .symbol = "Ac",
            .radius = static_cast<real_t>(1.86),
            .mass = static_cast<real_t>(227.0),
            .form_factors =
                {
                    static_cast<real_t>(15.914053),
                    static_cast<real_t>(0.080511),
                    static_cast<real_t>(32.535042),
                    static_cast<real_t>(0.770669),
                    static_cast<real_t>(21.553976),
                    static_cast<real_t>(4.352206),
                    static_cast<real_t>(11.433394),
                    static_cast<real_t>(21.381622),
                    static_cast<real_t>(3.939212),
                },
        },
        {
            .symbol = "Ag",
            .radius = static_cast<real_t>(1.28),
            .mass = static_cast<real_t>(107.87),
            .form_factors =
                {
                    static_cast<real_t>(6.073874),
                    static_cast<real_t>(0.055333),
                    static_cast<real_t>(17.155437),
                    static_cast<real_t>(7.896512),
                    static_cast<real_t>(4.173344),
                    static_cast<real_t>(28.443739),
                    static_cast<real_t>(0.852238),
                    static_cast<real_t>(110.37611),
                    static_cast<real_t>(0.756603),
                },
        },
        {
            .symbol = "Al",
            .radius = static_cast<real_t>(1.26),
            .mass = static_cast<real_t>(26.982),
            .form_factors =
                {
                    static_cast<real_t>(4.730796),
                    static_cast<real_t>(3.628931),
                    static_cast<real_t>(2.313951),
                    static_cast<real_t>(43.051167),
                    static_cast<real_t>(1.54198),
                    static_cast<real_t>(0.09596),
                    static_cast<real_t>(1.117564),
                    static_cast<real_t>(108.93239),
                    static_cast<real_t>(0.139509),
                },
        },
        {
            .symbol = "Am",
            .radius = static_cast<real_t>(1.66),
            .mass = static_cast<real_t>(243.0),
            .form_factors =
                {
                    static_cast<real_t>(33.435162),
                    static_cast<real_t>(0.612785),
                    static_cast<real_t>(23.657259),
                    static_cast<real_t>(3.792942),
                    static_cast<real_t>(15.576339),
                    static_cast<real_t>(16.195778),
                    static_cast<real_t>(3.027023),
                    static_cast<real_t>(117.757),
                    static_cast<real_t>(3.54116),
                },
        },
        {
            .symbol = "Ar",
            .radius = static_cast<real_t>(0.96),
            .mass = static_cast<real_t>(39.948),
            .form_factors =
                {
                    static_cast<real_t>(7.188004),
                    static_cast<real_t>(0.956221),
                    static_cast<real_t>(6.638454),
                    static_cast<real_t>(15.339877),
                    static_cast<real_t>(0.45418),
                    static_cast<real_t>(15.339862),
                    static_cast<real_t>(1.929593),
                    static_cast<real_t>(39.043823),
                    static_cast<real_t>(0.265954),
                },
        },
        {
            .symbol = "As",
            .radius = static_cast<real_t>(1.21),
            .mass = static_cast<real_t>(74.922),
            .form_factors =
                {
                    static_cast<real_t>(17.025642),
                    static_cast<real_t>(2.597739),
                    static_cast<real_t>(4.503441),
                    static_cast<real_t>(0.003012),
                    static_cast<real_t>(3.715904),
                    static_cast<real_t>(14.272119),
                    static_cast<real_t>(3.9372),
                    static_cast<real_t>(50.437996),
                    static_cast<real_t>(-2.984117),
                },
        },
        {
            .symbol = "At",
            .radius = static_cast<real_t>(1.47),
            .mass = static_cast<real_t>(210.0),
            .form_factors =
                {
                    static_cast<real_t>(16.011461),
                    static_cast<real_t>(0.092639),
                    static_cast<real_t>(32.615547),
                    static_cast<real_t>(0.904416),
                    static_cast<real_t>(8.113899),
                    static_cast<real_t>(26.543257),
                    static_cast<real_t>(2.884082),
                    static_cast<real_t>(68.372963),
                    static_cast<real_t>(3.995684),
                },
        },
        {
            .symbol = "Au",
            .radius = static_cast<real_t>(1.24),
            .mass = static_cast<real_t>(196.97),
            .form_factors =
                {
                    static_cast<real_t>(16.77739),
                    static_cast<real_t>(0.122737),
                    static_cast<real_t>(19.317156),
                    static_cast<real_t>(8.62157),
                    static_cast<real_t>(32.979683),
                    static_cast<real_t>(1.256902),
                    static_cast<real_t>(5.595453),
                    static_cast<real_t>(38.00882),
                    static_cast<real_t>(-6.279078),
                },
        },
        {
            .symbol = "B",
            .radius = static_cast<real_t>(0.85),
            .mass = static_cast<real_t>(10.81),
            .form_factors =
                {
                    static_cast<real_t>(2.085185),
                    static_cast<real_t>(23.494068),
                    static_cast<real_t>(1.06458),
                    static_cast<real_t>(1.137894),
                    static_cast<real_t>(1.062788),
                    static_cast<real_t>(61.238976),
                    static_cast<real_t>(0.140515),
                    static_cast<real_t>(0.114886),
                    static_cast<real_t>(0.003823),
                },
        },
        {
            .symbol = "Ba",
            .radius = static_cast<real_t>(1.96),
            .mass = static_cast<real_t>(137.33),
            .form_factors =
                {
                    static_cast<real_t>(19.747343),
                    static_cast<real_t>(3.481823),
                    static_cast<real_t>(17.368477),
                    static_cast<real_t>(0.371224),
                    static_cast<real_t>(10.465718),
                    static_cast<real_t>(21.226641),
                    static_cast<real_t>(2.592602),
                    static_cast<real_t>(173.83427),
                    static_cast<real_t>(-5.183497),
                },
        },
        {
            .symbol = "Be",
            .radius = static_cast<real_t>(1.02),
            .mass = static_cast<real_t>(9.0122),
            .form_factors =
                {
                    static_cast<real_t>(1.533712),
                    static_cast<real_t>(42.662079),
                    static_cast<real_t>(0.638283),
                    static_cast<real_t>(0.59542),
                    static_cast<real_t>(0.601052),
                    static_cast<real_t>(99.106499),
                    static_cast<real_t>(0.106139),
                    static_cast<real_t>(0.15134),
                    static_cast<real_t>(0.002511),
                },
        },
        {
            .symbol = "Bi",
            .radius = static_cast<real_t>(1.51),
            .mass = static_cast<real_t>(208.98),
            .form_factors =
                {
                    static_cast<real_t>(16.282274),
                    static_cast<real_t>(0.10118),
                    static_cast<real_t>(32.725136),
                    static_cast<real_t>(1.002287),
                    static_cast<real_t>(6.678302),
                    static_cast<real_t>(25.714146),
                    static_cast<real_t>(2.69475),
                    static_cast<real_t>(77.057549),
                    static_cast<real_t>(4.040914),
                },
        },
        {
            .symbol = "Bk",
            .radius = static_cast<real_t>(1.68),
            .mass = static_cast<real_t>(247.0),
            .form_factors =
                {
                    static_cast<real_t>(15.889072),
                    static_cast<real_t>(0.055503),
                    static_cast<real_t>(33.625286),
                    static_cast<real_t>(0.569571),
                    static_cast<real_t>(24.710381),
                    static_cast<real_t>(3.615472),
                    static_cast<real_t>(3.707139),
                    static_cast<real_t>(97.694786),
                    static_cast<real_t>(3.213169),
                },
        },
        {
            .symbol = "Br",
            .radius = static_cast<real_t>(1.14),
            .mass = static_cast<real_t>(79.904),
            .form_factors =
                {
                    static_cast<real_t>(17.55057),
                    static_cast<real_t>(2.119226),
                    static_cast<real_t>(5.411882),
                    static_cast<real_t>(16.557184),
                    static_cast<real_t>(3.93718),
                    static_cast<real_t>(0.002481),
                    static_cast<real_t>(3.880645),
                    static_cast<real_t>(42.164009),
                    static_cast<real_t>(-2.492088),
                },
        },
        {
            .symbol = "C",
            .radius = static_cast<real_t>(0.75),
            .mass = static_cast<real_t>(12.011),
            .form_factors =
                {
                    static_cast<real_t>(2.657506),
                    static_cast<real_t>(14.780758),
                    static_cast<real_t>(1.078079),
                    static_cast<real_t>(0.776775),
                    static_cast<real_t>(1.490909),
                    static_cast<real_t>(42.086842),
                    static_cast<real_t>(-4.24107),
                    static_cast<real_t>(-0.000294),
                    static_cast<real_t>(4.297983),
                },
        },
        {
            .symbol = "Ca",
            .radius = static_cast<real_t>(1.71),
            .mass = static_cast<real_t>(40.078),
            .form_factors =
                {
                    static_cast<real_t>(8.593655),
                    static_cast<real_t>(10.460644),
                    static_cast<real_t>(1.477324),
                    static_cast<real_t>(0.041891),
                    static_cast<real_t>(1.436254),
                    static_cast<real_t>(81.390381),
                    static_cast<real_t>(1.182839),
                    static_cast<real_t>(169.84784),
                    static_cast<real_t>(0.196255),
                },
        },
        {
            .symbol = "Cd",
            .radius = static_cast<real_t>(1.36),
            .mass = static_cast<real_t>(112.41),
            .form_factors =
                {
                    static_cast<real_t>(6.080986),
                    static_cast<real_t>(0.04899),
                    static_cast<real_t>(18.019468),
                    static_cast<real_t>(7.273646),
                    static_cast<real_t>(4.018197),
                    static_cast<real_t>(29.119284),
                    static_cast<real_t>(1.30351),
                    static_cast<real_t>(95.831207),
                    static_cast<real_t>(0.603504),
                },
        },
        {
            .symbol = "Ce",
            .radius = static_cast<real_t>(1.63),
            .mass = static_cast<real_t>(140.12),
            .form_factors =
                {
                    static_cast<real_t>(17.355122),
                    static_cast<real_t>(0.328369),
                    static_cast<real_t>(43.988499),
                    static_cast<real_t>(0.002047),
                    static_cast<real_t>(20.54665),
                    static_cast<real_t>(3.088196),
                    static_cast<real_t>(3.13067),
                    static_cast<real_t>(134.90765),
                    static_cast<real_t>(-38.386017),
                },
        },
        {
            .symbol = "Cf",
            .radius = static_cast<real_t>(1.68),
            .mass = static_cast<real_t>(251.0),
            .form_factors =
                {
                    static_cast<real_t>(33.794075),
                    static_cast<real_t>(0.550447),
                    static_cast<real_t>(25.467693),
                    static_cast<real_t>(3.581973),
                    static_cast<real_t>(16.048487),
                    static_cast<real_t>(14.357388),
                    static_cast<real_t>(3.657525),
                    static_cast<real_t>(96.064972),
                    static_cast<real_t>(3.005326),
                },
        },
        {
            .symbol = "Cl",
            .radius = static_cast<real_t>(0.99),
            .mass = static_cast<real_t>(35.45),
            .form_factors =
                {
                    static_cast<real_t>(1.446071),
                    static_cast<real_t>(0.052357),
                    static_cast<real_t>(6.870609),
                    static_cast<real_t>(1.193165),
                    static_cast<real_t>(6.151801),
                    static_cast<real_t>(18.343416),
                    static_cast<real_t>(1.750347),
                    static_cast<real_t>(46.398396),
                    static_cast<real_t>(0.146773),
                },
        },
        {
            .symbol = "Cm",
            .radius = static_cast<real_t>(1.66),
            .mass = static_cast<real_t>(247.0),
            .form_factors =
                {
                    static_cast<real_t>(15.804837),
                    static_cast<real_t>(0.058619),
                    static_cast<real_t>(33.480801),
                    static_cast<real_t>(0.59016),
                    static_cast<real_t>(24.150198),
                    static_cast<real_t>(3.67472),
                    static_cast<real_t>(3.655563),
                    static_cast<real_t>(100.73619),
                    static_cast<real_t>(3.39084),
                },
        },
        {
            .symbol = "Co",
            .radius = static_cast<real_t>(1.11),
            .mass = static_cast<real_t>(58.933),
            .form_factors =
                {
                    static_cast<real_t>(12.91451),
                    static_cast<real_t>(4.507138),
                    static_cast<real_t>(2.481908),
                    static_cast<real_t>(0.009126),
                    static_cast<real_t>(3.466894),
                    static_cast<real_t>(16.438129),
                    static_cast<real_t>(2.106351),
                    static_cast<real_t>(76.98732),
                    static_cast<real_t>(-0.936572),
                },
        },
        {
            .symbol = "Cr",
            .radius = static_cast<real_t>(1.22),
            .mass = static_cast<real_t>(51.996),
            .form_factors =
                {
                    static_cast<real_t>(11.007069),
                    static_cast<real_t>(6.366281),
                    static_cast<real_t>(1.555477),
                    static_cast<real_t>(0.023987),
                    static_cast<real_t>(2.985293),
                    static_cast<real_t>(23.244839),
                    static_cast<real_t>(1.347855),
                    static_cast<real_t>(105.7745),
                    static_cast<real_t>(0.06551),
                },
        },
        {
            .symbol = "Cs",
            .radius = static_cast<real_t>(2.32),
            .mass = static_cast<real_t>(132.91),
            .form_factors =
                {
                    static_cast<real_t>(17.418674),
                    static_cast<real_t>(0.399828),
                    static_cast<real_t>(8.314444),
                    static_cast<real_t>(0.016872),
                    static_cast<real_t>(10.323193),
                    static_cast<real_t>(25.605827),
                    static_cast<real_t>(1.383834),
                    static_cast<real_t>(233.33968),
                    static_cast<real_t>(-2.322802),
                },
        },
        {
            .symbol = "Cu",
            .radius = static_cast<real_t>(1.12),
            .mass = static_cast<real_t>(63.546),
            .form_factors =
                {
                    static_cast<real_t>(14.014192),
                    static_cast<real_t>(3.73828),
                    static_cast<real_t>(4.784577),
                    static_cast<real_t>(0.003744),
                    static_cast<real_t>(5.056806),
                    static_cast<real_t>(13.034982),
                    static_cast<real_t>(1.457971),
                    static_cast<real_t>(72.554794),
                    static_cast<real_t>(-3.254477),
                },
        },
        {
            .symbol = "Dy",
            .radius = static_cast<real_t>(1.67),
            .mass = static_cast<real_t>(162.50),
            .form_factors =
                {
                    static_cast<real_t>(26.671785),
                    static_cast<real_t>(2.282593),
                    static_cast<real_t>(88.687576),
                    static_cast<real_t>(0.000665),
                    static_cast<real_t>(14.065445),
                    static_cast<real_t>(12.92023),
                    static_cast<real_t>(2.768497),
                    static_cast<real_t>(121.93719),
                    static_cast<real_t>(-83.279831),
                },
        },
        {
            .symbol = "Er",
            .radius = static_cast<real_t>(1.65),
            .mass = static_cast<real_t>(167.26),
            .form_factors =
                {
                    static_cast<real_t>(28.174887),
                    static_cast<real_t>(2.120995),
                    static_cast<real_t>(82.493271),
                    static_cast<real_t>(0.00064),
                    static_cast<real_t>(14.624002),
                    static_cast<real_t>(11.915256),
                    static_cast<real_t>(2.802756),
                    static_cast<real_t>(114.52994),
                    static_cast<real_t>(-77.135223),
                },
        },
        {
            .symbol = "Eu",
            .radius = static_cast<real_t>(1.68),
            .mass = static_cast<real_t>(151.96),
            .form_factors =
                {
                    static_cast<real_t>(17.186195),
                    static_cast<real_t>(0.261678),
                    static_cast<real_t>(37.156837),
                    static_cast<real_t>(0.001995),
                    static_cast<real_t>(13.103387),
                    static_cast<real_t>(14.78736),
                    static_cast<real_t>(2.707246),
                    static_cast<real_t>(134.8163),
                    static_cast<real_t>(-31.586687),
                },
        },
        {
            .symbol = "F",
            .radius = static_cast<real_t>(0.64),
            .mass = static_cast<real_t>(18.998),
            .form_factors =
                {
                    static_cast<real_t>(3.511943),
                    static_cast<real_t>(10.687859),
                    static_cast<real_t>(2.772244),
                    static_cast<real_t>(4.380466),
                    static_cast<real_t>(0.678385),
                    static_cast<real_t>(0.093982),
                    static_cast<real_t>(0.915159),
                    static_cast<real_t>(27.255203),
                    static_cast<real_t>(0.032557),
                },
        },
        {
            .symbol = "Fe",
            .radius = static_cast<real_t>(1.16),
            .mass = static_cast<real_t>(55.845),
            .form_factors =
                {
                    static_cast<real_t>(12.311098),
                    static_cast<real_t>(5.009415),
                    static_cast<real_t>(1.876623),
                    static_cast<real_t>(0.014461),
                    static_cast<real_t>(3.066177),
                    static_cast<real_t>(18.74304),
                    static_cast<real_t>(2.070451),
                    static_cast<real_t>(82.767876),
                    static_cast<real_t>(-0.304931),
                },
        },
        {
            .symbol = "Fr",
            .radius = static_cast<real_t>(2.23),
            .mass = static_cast<real_t>(223.0),
            .form_factors =
                {
                    static_cast<real_t>(16.007385),
                    static_cast<real_t>(0.087031),
                    static_cast<real_t>(32.66383),
                    static_cast<real_t>(0.840187),
                    static_cast<real_t>(21.594351),
                    static_cast<real_t>(4.954467),
                    static_cast<real_t>(1.598497),
                    static_cast<real_t>(199.8058),
                    static_cast<real_t>(4.003472),
                },
        },
        {
            .symbol = "Ga",
            .radius = static_cast<real_t>(1.24),
            .mass = static_cast<real_t>(69.723),
            .form_factors =
                {
                    static_cast<real_t>(15.758946),
                    static_cast<real_t>(3.121754),
                    static_cast<real_t>(6.841123),
                    static_cast<real_t>(0.226057),
                    static_cast<real_t>(4.121016),
                    static_cast<real_t>(12.482196),
                    static_cast<real_t>(2.714681),
                    static_cast<real_t>(66.203621),
                    static_cast<real_t>(-0.847395),
                },
        },
        {
            .symbol = "Gd",
            .radius = static_cast<real_t>(1.69),
            .mass = static_cast<real_t>(157.25),
            .form_factors =
                {
                    static_cast<real_t>(24.898117),
                    static_cast<real_t>(2.435028),
                    static_cast<real_t>(17.104952),
                    static_cast<real_t>(0.246961),
                    static_cast<real_t>(13.222581),
                    static_cast<real_t>(13.996325),
                    static_cast<real_t>(3.266152),
                    static_cast<real_t>(110.86309),
                    static_cast<real_t>(-43.505684),
                },
        },
        {
            .symbol = "Ge",
            .radius = static_cast<real_t>(1.21),
            .mass = static_cast<real_t>(72.630),
            .form_factors =
                {
                    static_cast<real_t>(16.540613),
                    static_cast<real_t>(2.866618),
                    static_cast<real_t>(1.5679),
                    static_cast<real_t>(0.012198),
                    static_cast<real_t>(3.727829),
                    static_cast<real_t>(13.432163),
                    static_cast<real_t>(3.345098),
                    static_cast<real_t>(58.866047),
                    static_cast<real_t>(0.018726),
                },
        },
        {
            .symbol = "H",
            .radius = static_cast<real_t>(0.32),
            .mass = static_cast<real_t>(1.008),
            .form_factors =
                {
                    static_cast<real_t>(0.413048),
                    static_cast<real_t>(15.569946),
                    static_cast<real_t>(0.294953),
                    static_cast<real_t>(32.398468),
                    static_cast<real_t>(0.187491),
                    static_cast<real_t>(5.711404),
                    static_cast<real_t>(0.080701),
                    static_cast<real_t>(61.889874),
                    static_cast<real_t>(4.9e-05),
                },
        },
        {
            .symbol = "He",
            .radius = static_cast<real_t>(0.46),
            .mass = static_cast<real_t>(4.0026),
            .form_factors =
                {
                    static_cast<real_t>(0.732354),
                    static_cast<real_t>(11.553918),
                    static_cast<real_t>(0.753896),
                    static_cast<real_t>(4.595831),
                    static_cast<real_t>(0.283819),
                    static_cast<real_t>(1.546299),
                    static_cast<real_t>(0.190003),
                    static_cast<real_t>(26.463964),
                    static_cast<real_t>(0.000487),
                },
        },
        {
            .symbol = "Hf",
            .radius = static_cast<real_t>(1.52),
            .mass = static_cast<real_t>(178.49),
            .form_factors =
                {
                    static_cast<real_t>(30.617033),
                    static_cast<real_t>(1.795613),
                    static_cast<real_t>(15.145351),
                    static_cast<real_t>(9.934469),
                    static_cast<real_t>(54.933548),
                    static_cast<real_t>(0.000739),
                    static_cast<real_t>(4.096253),
                    static_cast<real_t>(76.189705),
                    static_cast<real_t>(-49.719837),
                },
        },
        {
            .symbol = "Hg",
            .radius = static_cast<real_t>(1.33),
            .mass = static_cast<real_t>(200.59),
            .form_factors =
                {
                    static_cast<real_t>(16.83989),
                    static_cast<real_t>(0.115905),
                    static_cast<real_t>(20.023823),
                    static_cast<real_t>(8.256927),
                    static_cast<real_t>(28.428564),
                    static_cast<real_t>(1.19525),
                    static_cast<real_t>(5.881564),
                    static_cast<real_t>(39.247227),
                    static_cast<real_t>(4.076478),
                },
        },
        {
            .symbol = "Ho",
            .radius = static_cast<real_t>(1.66),
            .mass = static_cast<real_t>(164.93),
            .form_factors =
                {
                    static_cast<real_t>(27.15019),
                    static_cast<real_t>(2.16966),
                    static_cast<real_t>(16.999819),
                    static_cast<real_t>(0.215414),
                    static_cast<real_t>(14.059334),
                    static_cast<real_t>(12.213148),
                    static_cast<real_t>(3.386979),
                    static_cast<real_t>(100.50678),
                    static_cast<real_t>(-41.165253),
                },
        },
        {
            .symbol = "I",
            .radius = static_cast<real_t>(1.33),
            .mass = static_cast<real_t>(126.90),
            .form_factors =
                {
                    static_cast<real_t>(19.884502),
                    static_cast<real_t>(4.628591),
                    static_cast<real_t>(6.736593),
                    static_cast<real_t>(0.027754),
                    static_cast<real_t>(8.110516),
                    static_cast<real_t>(31.849096),
                    static_cast<real_t>(1.170953),
                    static_cast<real_t>(84.406387),
                    static_cast<real_t>(-0.448811),
                },
        },
        {
            .symbol = "In",
            .radius = static_cast<real_t>(1.42),
            .mass = static_cast<real_t>(114.82),
            .form_factors =
                {
                    static_cast<real_t>(6.196477),
                    static_cast<real_t>(0.042072),
                    static_cast<real_t>(18.816183),
                    static_cast<real_t>(6.695665),
                    static_cast<real_t>(4.050479),
                    static_cast<real_t>(31.00979),
                    static_cast<real_t>(1.638929),
                    static_cast<real_t>(103.28435),
                    static_cast<real_t>(0.333097),
                },
        },
        {
            .symbol = "Ir",
            .radius = static_cast<real_t>(1.22),
            .mass = static_cast<real_t>(192.22),
            .form_factors =
                {
                    static_cast<real_t>(32.004436),
                    static_cast<real_t>(1.353767),
                    static_cast<real_t>(1.975454),
                    static_cast<real_t>(81.014175),
                    static_cast<real_t>(17.070105),
                    static_cast<real_t>(0.128093),
                    static_cast<real_t>(15.939454),
                    static_cast<real_t>(7.661196),
                    static_cast<real_t>(4.018893),
                },
        },
        {
            .symbol = "K",
            .radius = static_cast<real_t>(1.96),
            .mass = static_cast<real_t>(39.098),
            .form_factors =
                {
                    static_cast<real_t>(8.163991),
                    static_cast<real_t>(12.816323),
                    static_cast<real_t>(7.146945),
                    static_cast<real_t>(0.808945),
                    static_cast<real_t>(1.07014),
                    static_cast<real_t>(210.32701),
                    static_cast<real_t>(0.877316),
                    static_cast<real_t>(39.597652),
                    static_cast<real_t>(0.253614),
                },
        },
        {
            .symbol = "Kr",
            .radius = static_cast<real_t>(1.17),
            .mass = static_cast<real_t>(83.798),
            .form_factors =
                {
                    static_cast<real_t>(17.655279),
                    static_cast<real_t>(1.908231),
                    static_cast<real_t>(6.848105),
                    static_cast<real_t>(16.606236),
                    static_cast<real_t>(4.171004),
                    static_cast<real_t>(0.001598),
                    static_cast<real_t>(3.44676),
                    static_cast<real_t>(39.917473),
                    static_cast<real_t>(-2.810592),
                },
        },
        {
            .symbol = "La",
            .radius = static_cast<real_t>(1.80),
            .mass = static_cast<real_t>(138.91),
            .form_factors =
                {
                    static_cast<real_t>(19.966019),
                    static_cast<real_t>(3.197408),
                    static_cast<real_t>(27.329655),
                    static_cast<real_t>(0.003446),
                    static_cast<real_t>(11.018425),
                    static_cast<real_t>(19.955492),
                    static_cast<real_t>(3.086696),
                    static_cast<real_t>(141.38197),
                    static_cast<real_t>(-21.745489),
                },
        },
        {
            .symbol = "Li",
            .radius = static_cast<real_t>(1.33),
            .mass = static_cast<real_t>(6.94),
            .form_factors =
                {
                    static_cast<real_t>(0.974637),
                    static_cast<real_t>(4.334946),
                    static_cast<real_t>(0.158472),
                    static_cast<real_t>(0.342451),
                    static_cast<real_t>(0.811855),
                    static_cast<real_t>(97.102966),
                    static_cast<real_t>(0.262416),
                    static_cast<real_t>(201.36383),
                    static_cast<real_t>(0.002542),
                },
        },
        {
            .symbol = "Lu",
            .radius = static_cast<real_t>(1.62),
            .mass = static_cast<real_t>(174.97),
            .form_factors =
                {
                    static_cast<real_t>(30.122866),
                    static_cast<real_t>(1.88309),
                    static_cast<real_t>(15.099346),
                    static_cast<real_t>(10.342764),
                    static_cast<real_t>(56.314899),
                    static_cast<real_t>(0.00078),
                    static_cast<real_t>(3.54098),
                    static_cast<real_t>(89.55925),
                    static_cast<real_t>(-51.049416),
                },
        },
        {
            .symbol = "Mg",
            .radius = static_cast<real_t>(1.39),
            .mass = static_cast<real_t>(24.305),
            .form_factors =
                {
                    static_cast<real_t>(4.708971),
                    static_cast<real_t>(4.875207),
                    static_cast<real_t>(1.194814),
                    static_cast<real_t>(108.50608),
                    static_cast<real_t>(1.558157),
                    static_cast<real_t>(0.111516),
                    static_cast<real_t>(1.170413),
                    static_cast<real_t>(48.292408),
                    static_cast<real_t>(0.126842),
                },
        },
        {
            .symbol = "Mn",
            .radius = static_cast<real_t>(1.19),
            .mass = static_cast<real_t>(54.938),
            .form_factors =
                {
                    static_cast<real_t>(11.709542),
                    static_cast<real_t>(5.59712),
                    static_cast<real_t>(1.733414),
                    static_cast<real_t>(0.0178),
                    static_cast<real_t>(2.673141),
                    static_cast<real_t>(21.78842),
                    static_cast<real_t>(2.023368),
                    static_cast<real_t>(89.517914),
                    static_cast<real_t>(-0.147293),
                },
        },
        {
            .symbol = "Mo",
            .radius = static_cast<real_t>(1.38),
            .mass = static_cast<real_t>(95.95),
            .form_factors =
                {
                    static_cast<real_t>(6.236218),
                    static_cast<real_t>(0.09078),
                    static_cast<real_t>(17.987711),
                    static_cast<real_t>(1.10831),
                    static_cast<real_t>(12.973127),
                    static_cast<real_t>(11.46872),
                    static_cast<real_t>(3.451426),
                    static_cast<real_t>(66.684151),
                    static_cast<real_t>(1.10877),
                },
        },
        {
            .symbol = "N",
            .radius = static_cast<real_t>(0.71),
            .mass = static_cast<real_t>(14.007),
            .form_factors =
                {
                    static_cast<real_t>(11.89378),
                    static_cast<real_t>(0.000158),
                    static_cast<real_t>(3.277479),
                    static_cast<real_t>(10.232723),
                    static_cast<real_t>(1.858092),
                    static_cast<real_t>(30.34469),
                    static_cast<real_t>(0.858927),
                    static_cast<real_t>(0.656065),
                    static_cast<real_t>(-11.804902),
                },
        },
        {
            .symbol = "Na",
            .radius = static_cast<real_t>(1.55),
            .mass = static_cast<real_t>(22.990),
            .form_factors =
                {
                    static_cast<real_t>(4.910127),
                    static_cast<real_t>(3.281434),
                    static_cast<real_t>(3.081783),
                    static_cast<real_t>(9.119178),
                    static_cast<real_t>(1.262067),
                    static_cast<real_t>(0.102763),
                    static_cast<real_t>(1.098938),
                    static_cast<real_t>(132.01395),
                    static_cast<real_t>(0.079712),
                },
        },
        {
            .symbol = "Nb",
            .radius = static_cast<real_t>(1.47),
            .mass = static_cast<real_t>(92.906),
            .form_factors =
                {
                    static_cast<real_t>(17.958399),
                    static_cast<real_t>(1.21159),
                    static_cast<real_t>(12.063054),
                    static_cast<real_t>(12.246687),
                    static_cast<real_t>(5.007015),
                    static_cast<real_t>(0.098615),
                    static_cast<real_t>(3.287667),
                    static_cast<real_t>(75.011948),
                    static_cast<real_t>(1.123452),
                },
        },
        {
            .symbol = "Nd",
            .radius = static_cast<real_t>(1.74),
            .mass = static_cast<real_t>(144.24),
            .form_factors =
                {
                    static_cast<real_t>(17.331244),
                    static_cast<real_t>(0.300269),
                    static_cast<real_t>(62.783924),
                    static_cast<real_t>(0.00132),
                    static_cast<real_t>(12.160097),
                    static_cast<real_t>(17.026001),
                    static_cast<real_t>(2.663483),
                    static_cast<real_t>(148.74899),
                    static_cast<real_t>(-57.189842),
                },
        },
        {
            .symbol = "Ne",
            .radius = static_cast<real_t>(0.67),
            .mass = static_cast<real_t>(20.180),
            .form_factors =
                {
                    static_cast<real_t>(4.183749),
                    static_cast<real_t>(8.175457),
                    static_cast<real_t>(2.905726),
                    static_cast<real_t>(3.252536),
                    static_cast<real_t>(0.520513),
                    static_cast<real_t>(0.063295),
                    static_cast<real_t>(1.135641),
                    static_cast<real_t>(21.81391),
                    static_cast<real_t>(0.025576),
                },
        },
        {
            .symbol = "Ni",
            .radius = static_cast<real_t>(1.10),
            .mass = static_cast<real_t>(58.693),
            .form_factors =
                {
                    static_cast<real_t>(13.521865),
                    static_cast<real_t>(4.077277),
                    static_cast<real_t>(6.947285),
                    static_cast<real_t>(0.286763),
                    static_cast<real_t>(3.866028),
                    static_cast<real_t>(14.622634),
                    static_cast<real_t>(2.1359),
                    static_cast<real_t>(71.96608),
                    static_cast<real_t>(-2.762697),
                },
        },
        {
            .symbol = "Np",
            .radius = static_cast<real_t>(1.71),
            .mass = static_cast<real_t>(237.0),
            .form_factors =
                {
                    static_cast<real_t>(32.999901),
                    static_cast<real_t>(0.657086),
                    static_cast<real_t>(22.638077),
                    static_cast<real_t>(3.854918),
                    static_cast<real_t>(14.219973),
                    static_cast<real_t>(17.435474),
                    static_cast<real_t>(3.67295),
                    static_cast<real_t>(109.46448),
                    static_cast<real_t>(3.769391),
                },
        },
        {
            .symbol = "O",
            .radius = static_cast<real_t>(0.63),
            .mass = static_cast<real_t>(15.999),
            .form_factors =
                {
                    static_cast<real_t>(2.960427),
                    static_cast<real_t>(14.182259),
                    static_cast<real_t>(2.508818),
                    static_cast<real_t>(5.936858),
                    static_cast<real_t>(0.637853),
                    static_cast<real_t>(0.112726),
                    static_cast<real_t>(0.722838),
                    static_cast<real_t>(34.958481),
                    static_cast<real_t>(0.027014),
                },
        },
        {
            .symbol = "Os",
            .radius = static_cast<real_t>(1.29),
            .mass = static_cast<real_t>(190.23),
            .form_factors =
                {
                    static_cast<real_t>(32.210297),
                    static_cast<real_t>(1.473531),
                    static_cast<real_t>(16.67844),
                    static_cast<real_t>(9.049695),
                    static_cast<real_t>(48.559906),
                    static_cast<real_t>(0.000519),
                    static_cast<real_t>(5.455839),
                    static_cast<real_t>(50.210201),
                    static_cast<real_t>(-43.677956),
                },
        },
        {
            .symbol = "P",
            .radius = static_cast<real_t>(1.11),
            .mass = static_cast<real_t>(30.974),
            .form_factors =
                {
                    static_cast<real_t>(1.950541),
                    static_cast<real_t>(0.908139),
                    static_cast<real_t>(4.14693),
                    static_cast<real_t>(27.044952),
                    static_cast<real_t>(1.49456),
                    static_cast<real_t>(0.07128),
                    static_cast<real_t>(1.522042),
                    static_cast<real_t>(67.520187),
                    static_cast<real_t>(0.155233),
                },
        },
        {
            .symbol = "Pa",
            .radius = static_cast<real_t>(1.69),
            .mass = static_cast<real_t>(231.04),
            .form_factors =
                {
                    static_cast<real_t>(32.740208),
                    static_cast<real_t>(0.709545),
                    static_cast<real_t>(21.973675),
                    static_cast<real_t>(4.050881),
                    static_cast<real_t>(12.957398),
                    static_cast<real_t>(19.231543),
                    static_cast<real_t>(3.683832),
                    static_cast<real_t>(117.255),
                    static_cast<real_t>(3.886066),
                },
        },
        {
            .symbol = "Pb",
            .radius = static_cast<real_t>(1.44),
            .mass = static_cast<real_t>(207.2),
            .form_factors =
                {
                    static_cast<real_t>(16.419567),
                    static_cast<real_t>(0.105499),
                    static_cast<real_t>(32.73859),
                    static_cast<real_t>(1.055049),
                    static_cast<real_t>(6.530247),
                    static_cast<real_t>(25.02589),
                    static_cast<real_t>(2.342742),
                    static_cast<real_t>(80.906593),
                    static_cast<real_t>(4.049824),
                },
        },
        {
            .symbol = "Pd",
            .radius = static_cast<real_t>(1.20),
            .mass = static_cast<real_t>(106.42),
            .form_factors =
                {
                    static_cast<real_t>(6.121511),
                    static_cast<real_t>(0.062549),
                    static_cast<real_t>(4.784063),
                    static_cast<real_t>(0.784031),
                    static_cast<real_t>(16.631683),
                    static_cast<real_t>(8.751391),
                    static_cast<real_t>(4.318258),
                    static_cast<real_t>(34.489983),
                    static_cast<real_t>(0.883099),
                },
        },
        {
            .symbol = "Pm",
            .radius = static_cast<real_t>(1.73),
            .mass = static_cast<real_t>(145.0),
            .form_factors =
                {
                    static_cast<real_t>(17.286388),
                    static_cast<real_t>(0.28662),
                    static_cast<real_t>(51.560162),
                    static_cast<real_t>(0.00155),
                    static_cast<real_t>(12.478557),
                    static_cast<real_t>(16.223755),
                    static_cast<real_t>(2.675515),
                    static_cast<real_t>(143.98451),
                    static_cast<real_t>(-45.973682),
                },
        },
        {
            .symbol = "Po",
            .radius = static_cast<real_t>(1.45),
            .mass = static_cast<real_t>(209.0),
            .form_factors =
                {
                    static_cast<real_t>(16.289164),
                    static_cast<real_t>(0.098121),
                    static_cast<real_t>(32.807171),
                    static_cast<real_t>(0.966265),
                    static_cast<real_t>(21.095163),
                    static_cast<real_t>(6.046622),
                    static_cast<real_t>(2.505901),
                    static_cast<real_t>(76.598068),
                    static_cast<real_t>(4.046556),
                },
        },
        {
            .symbol = "Pr",
            .radius = static_cast<real_t>(1.76),
            .mass = static_cast<real_t>(140.91),
            .form_factors =
                {
                    static_cast<real_t>(21.551311),
                    static_cast<real_t>(2.995675),
                    static_cast<real_t>(17.16173),
                    static_cast<real_t>(0.312491),
                    static_cast<real_t>(11.903859),
                    static_cast<real_t>(17.716705),
                    static_cast<real_t>(2.679103),
                    static_cast<real_t>(152.19282),
                    static_cast<real_t>(-3.871068),
                },
        },
        {
            .symbol = "Pt",
            .radius = static_cast<real_t>(1.23),
            .mass = static_cast<real_t>(195.08),
            .form_factors =
                {
                    static_cast<real_t>(31.273891),
                    static_cast<real_t>(1.316992),
                    static_cast<real_t>(18.44544),
                    static_cast<real_t>(8.797154),
                    static_cast<real_t>(17.063745),
                    static_cast<real_t>(0.124741),
                    static_cast<real_t>(5.555933),
                    static_cast<real_t>(40.177994),
                    static_cast<real_t>(4.050394),
                },
        },
        {
            .symbol = "Pu",
            .radius = static_cast<real_t>(1.72),
            .mass = static_cast<real_t>(244.0),
            .form_factors =
                {
                    static_cast<real_t>(33.281178),
                    static_cast<real_t>(0.634999),
                    static_cast<real_t>(23.148544),
                    static_cast<real_t>(3.856168),
                    static_cast<real_t>(15.153755),
                    static_cast<real_t>(16.849735),
                    static_cast<real_t>(3.031492),
                    static_cast<real_t>(121.29204),
                    static_cast<real_t>(3.6642),
                },
        },
        {
            .symbol = "Ra",
            .radius = static_cast<real_t>(2.01),
            .mass = static_cast<real_t>(226.0),
            .form_factors =
                {
                    static_cast<real_t>(32.56369),
                    static_cast<real_t>(0.80198),
                    static_cast<real_t>(21.396671),
                    static_cast<real_t>(4.590666),
                    static_cast<real_t>(11.298093),
                    static_cast<real_t>(22.758972),
                    static_cast<real_t>(2.834688),
                    static_cast<real_t>(160.40439),
                    static_cast<real_t>(3.981773),
                },
        },
        {
            .symbol = "Rb",
            .radius = static_cast<real_t>(2.10),
            .mass = static_cast<real_t>(85.468),
            .form_factors =
                {
                    static_cast<real_t>(8.123134),
                    static_cast<real_t>(15.142385),
                    static_cast<real_t>(2.138042),
                    static_cast<real_t>(33.542667),
                    static_cast<real_t>(6.761702),
                    static_cast<real_t>(0.129372),
                    static_cast<real_t>(1.156051),
                    static_cast<real_t>(224.13251),
                    static_cast<real_t>(1.139548),
                },
        },
        {
            .symbol = "Re",
            .radius = static_cast<real_t>(1.31),
            .mass = static_cast<real_t>(186.21),
            .form_factors =
                {
                    static_cast<real_t>(31.888456),
                    static_cast<real_t>(1.549238),
                    static_cast<real_t>(16.117104),
                    static_cast<real_t>(9.233474),
                    static_cast<real_t>(42.390297),
                    static_cast<real_t>(0.000689),
                    static_cast<real_t>(5.211669),
                    static_cast<real_t>(54.516373),
                    static_cast<real_t>(-37.412682),
                },
        },
        {
            .symbol = "Rh",
            .radius = static_cast<real_t>(1.25),
            .mass = static_cast<real_t>(102.91),
            .form_factors =
                {
                    static_cast<real_t>(6.216648),
                    static_cast<real_t>(0.070789),
                    static_cast<real_t>(17.919739),
                    static_cast<real_t>(0.856121),
                    static_cast<real_t>(3.854252),
                    static_cast<real_t>(33.889484),
                    static_cast<real_t>(0.840326),
                    static_cast<real_t>(121.68669),
                    static_cast<real_t>(0.995452),
                },
        },
        {
            .symbol = "Rn",
            .radius = static_cast<real_t>(1.42),
            .mass = static_cast<real_t>(222.0),
            .form_factors =
                {
                    static_cast<real_t>(16.070229),
                    static_cast<real_t>(0.090437),
                    static_cast<real_t>(32.641106),
                    static_cast<real_t>(0.876409),
                    static_cast<real_t>(21.489658),
                    static_cast<real_t>(5.239687),
                    static_cast<real_t>(2.299218),
                    static_cast<real_t>(69.188477),
                    static_cast<real_t>(4.020977),
                },
        },
        {
            .symbol = "Ru",
            .radius = static_cast<real_t>(1.25),
            .mass = static_cast<real_t>(101.07),
            .form_factors =
                {
                    static_cast<real_t>(6.271624),
                    static_cast<real_t>(0.07704),
                    static_cast<real_t>(17.906738),
                    static_cast<real_t>(0.928222),
                    static_cast<real_t>(14.123269),
                    static_cast<real_t>(9.555345),
                    static_cast<real_t>(3.746008),
                    static_cast<real_t>(35.86068),
                    static_cast<real_t>(1.043992),
                },
        },
        {
            .symbol = "S",
            .radius = static_cast<real_t>(1.03),
            .mass = static_cast<real_t>(32.06),
            .form_factors =
                {
                    static_cast<real_t>(6.372157),
                    static_cast<real_t>(1.514347),
                    static_cast<real_t>(5.154568),
                    static_cast<real_t>(22.092527),
                    static_cast<real_t>(1.473732),
                    static_cast<real_t>(0.061373),
                    static_cast<real_t>(1.635073),
                    static_cast<real_t>(55.445175),
                    static_cast<real_t>(0.154722),
                },
        },
        {
            .symbol = "Sb",
            .radius = static_cast<real_t>(1.40),
            .mass = static_cast<real_t>(121.76),
            .form_factors =
                {
                    static_cast<real_t>(5.394956),
                    static_cast<real_t>(33.326523),
                    static_cast<real_t>(6.54957),
                    static_cast<real_t>(0.030974),
                    static_cast<real_t>(19.650681),
                    static_cast<real_t>(5.564929),
                    static_cast<real_t>(1.82782),
                    static_cast<real_t>(87.130966),
                    static_cast<real_t>(-0.290506),
                },
        },
        {
            .symbol = "Sc",
            .radius = static_cast<real_t>(1.48),
            .mass = static_cast<real_t>(44.956),
            .form_factors =
                {
                    static_cast<real_t>(1.476566),
                    static_cast<real_t>(53.131023),
                    static_cast<real_t>(1.487278),
                    static_cast<real_t>(0.035325),
                    static_cast<real_t>(1.600187),
                    static_cast<real_t>(137.31949),
                    static_cast<real_t>(9.177463),
                    static_cast<real_t>(9.098031),
                    static_cast<real_t>(0.157765),
                },
        },
        {
            .symbol = "Se",
            .radius = static_cast<real_t>(1.16),
            .mass = static_cast<real_t>(78.971),
            .form_factors =
                {
                    static_cast<real_t>(17.354071),
                    static_cast<real_t>(2.349787),
                    static_cast<real_t>(4.653248),
                    static_cast<real_t>(0.00255),
                    static_cast<real_t>(4.259489),
                    static_cast<real_t>(15.57946),
                    static_cast<real_t>(4.136455),
                    static_cast<real_t>(45.181202),
                    static_cast<real_t>(-3.160982),
                },
        },
        {
            .symbol = "Si",
            .radius = static_cast<real_t>(1.16),
            .mass = static_cast<real_t>(28.085),
            .form_factors =
                {
                    static_cast<real_t>(5.275329),
                    static_cast<real_t>(2.631338),
                    static_cast<real_t>(3.191038),
                    static_cast<real_t>(33.730728),
                    static_cast<real_t>(1.511514),
                    static_cast<real_t>(0.081119),
                    static_cast<real_t>(1.356849),
                    static_cast<real_t>(86.288643),
                    static_cast<real_t>(0.145073),
                },
        },
        {
            .symbol = "Sm",
            .radius = static_cast<real_t>(1.72),
            .mass = static_cast<real_t>(150.36),
            .form_factors =
                {
                    static_cast<real_t>(23.700363),
                    static_cast<real_t>(2.689539),
                    static_cast<real_t>(23.072214),
                    static_cast<real_t>(0.003491),
                    static_cast<real_t>(12.777782),
                    static_cast<real_t>(15.495437),
                    static_cast<real_t>(2.684217),
                    static_cast<real_t>(139.86247),
                    static_cast<real_t>(-17.452166),
                },
        },
        {
            .symbol = "Sn",
            .radius = static_cast<real_t>(1.40),
            .mass = static_cast<real_t>(118.71),
            .form_factors =
                {
                    static_cast<real_t>(19.325171),
                    static_cast<real_t>(6.118104),
                    static_cast<real_t>(6.281571),
                    static_cast<real_t>(0.036915),
                    static_cast<real_t>(4.498866),
                    static_cast<real_t>(32.529045),
                    static_cast<real_t>(1.856934),
                    static_cast<real_t>(95.037186),
                    static_cast<real_t>(0.119024),
                },
        },
        {
            .symbol = "Sr",
            .radius = static_cast<real_t>(1.85),
            .mass = static_cast<real_t>(87.62),
            .form_factors =
                {
                    static_cast<real_t>(17.730219),
                    static_cast<real_t>(1.56306),
                    static_cast<real_t>(9.795867),
                    static_cast<real_t>(14.310868),
                    static_cast<real_t>(6.099763),
                    static_cast<real_t>(0.120574),
                    static_cast<real_t>(2.620025),
                    static_cast<real_t>(135.77132),
                    static_cast<real_t>(1.140251),
                },
        },
        {
            .symbol = "Ta",
            .radius = static_cast<real_t>(1.46),
            .mass = static_cast<real_t>(180.95),
            .form_factors =
                {
                    static_cast<real_t>(31.066359),
                    static_cast<real_t>(1.708732),
                    static_cast<real_t>(15.341823),
                    static_cast<real_t>(9.618455),
                    static_cast<real_t>(49.278297),
                    static_cast<real_t>(0.00076),
                    static_cast<real_t>(4.577665),
                    static_cast<real_t>(66.346199),
                    static_cast<real_t>(-44.119026),
                },
        },
        {
            .symbol = "Tb",
            .radius = static_cast<real_t>(1.68),
            .mass = static_cast<real_t>(158.93),
            .form_factors =
                {
                    static_cast<real_t>(25.910013),
                    static_cast<real_t>(2.373912),
                    static_cast<real_t>(32.344139),
                    static_cast<real_t>(0.002034),
                    static_cast<real_t>(13.765117),
                    static_cast<real_t>(13.481969),
                    static_cast<real_t>(2.751404),
                    static_cast<real_t>(125.83651),
                    static_cast<real_t>(-26.851971),
                },
        },
        {
            .symbol = "Tc",
            .radius = static_cast<real_t>(1.28),
            .mass = static_cast<real_t>(98.0),
            .form_factors =
                {
                    static_cast<real_t>(17.840963),
                    static_cast<real_t>(1.005729),
                    static_cast<real_t>(3.428236),
                    static_cast<real_t>(41.901382),
                    static_cast<real_t>(1.373012),
                    static_cast<real_t>(119.32054),
                    static_cast<real_t>(12.947364),
                    static_cast<real_t>(9.781542),
                    static_cast<real_t>(1.074784),
                },
        },
        {
            .symbol = "Te",
            .radius = static_cast<real_t>(1.36),
            .mass = static_cast<real_t>(127.60),
            .form_factors =
                {
                    static_cast<real_t>(6.660302),
                    static_cast<real_t>(33.031654),
                    static_cast<real_t>(6.940756),
                    static_cast<real_t>(0.02575),
                    static_cast<real_t>(19.847015),
                    static_cast<real_t>(5.065547),
                    static_cast<real_t>(1.557175),
                    static_cast<real_t>(84.101616),
                    static_cast<real_t>(-0.806668),
                },
        },
        {
            .symbol = "Th",
            .radius = static_cast<real_t>(1.75),
            .mass = static_cast<real_t>(232.04),
            .form_factors =
                {
                    static_cast<real_t>(15.784024),
                    static_cast<real_t>(0.077067),
                    static_cast<real_t>(32.454899),
                    static_cast<real_t>(0.735137),
                    static_cast<real_t>(21.849222),
                    static_cast<real_t>(4.097976),
                    static_cast<real_t>(4.239077),
                    static_cast<real_t>(109.46411),
                    static_cast<real_t>(3.922533),
                },
        },
        {
            .symbol = "Ti",
            .radius = static_cast<real_t>(1.36),
            .mass = static_cast<real_t>(47.867),
            .form_factors =
                {
                    static_cast<real_t>(9.818524),
                    static_cast<real_t>(8.001879),
                    static_cast<real_t>(1.522646),
                    static_cast<real_t>(0.029763),
                    static_cast<real_t>(1.703101),
                    static_cast<real_t>(39.885422),
                    static_cast<real_t>(1.768774),
                    static_cast<real_t>(120.158),
                    static_cast<real_t>(0.102473),
                },
        },
        {
            .symbol = "Tl",
            .radius = static_cast<real_t>(1.44),
            .mass = static_cast<real_t>(204.38),
            .form_factors =
                {
                    static_cast<real_t>(16.630795),
                    static_cast<real_t>(0.110704),
                    static_cast<real_t>(19.386616),
                    static_cast<real_t>(7.181401),
                    static_cast<real_t>(32.808571),
                    static_cast<real_t>(1.11973),
                    static_cast<real_t>(1.747191),
                    static_cast<real_t>(90.660263),
                    static_cast<real_t>(4.066939),
                },
        },
        {
            .symbol = "Tm",
            .radius = static_cast<real_t>(1.64),
            .mass = static_cast<real_t>(168.93),
            .form_factors =
                {
                    static_cast<real_t>(28.925894),
                    static_cast<real_t>(2.046203),
                    static_cast<real_t>(76.173798),
                    static_cast<real_t>(0.000656),
                    static_cast<real_t>(14.904704),
                    static_cast<real_t>(11.465375),
                    static_cast<real_t>(2.814812),
                    static_cast<real_t>(111.41198),
                    static_cast<real_t>(-70.839813),
                },
        },
        {
            .symbol = "U",
            .radius = static_cast<real_t>(1.70),
            .mass = static_cast<real_t>(238.03),
            .form_factors =
                {
                    static_cast<real_t>(15.679275),
                    static_cast<real_t>(0.071206),
                    static_cast<real_t>(32.824306),
                    static_cast<real_t>(0.681177),
                    static_cast<real_t>(13.660459),
                    static_cast<real_t>(18.236156),
                    static_cast<real_t>(3.687261),
                    static_cast<real_t>(112.50004),
                    static_cast<real_t>(3.854444),
                },
        },
        {
            .symbol = "V",
            .radius = static_cast<real_t>(1.34),
            .mass = static_cast<real_t>(50.942),
            .form_factors =
                {
                    static_cast<real_t>(10.473575),
                    static_cast<real_t>(7.08194),
                    static_cast<real_t>(1.547881),
                    static_cast<real_t>(0.02604),
                    static_cast<real_t>(1.986381),
                    static_cast<real_t>(31.909672),
                    static_cast<real_t>(1.865616),
                    static_cast<real_t>(108.02284),
                    static_cast<real_t>(0.067744),
                },
        },
        {
            .symbol = "W",
            .radius = static_cast<real_t>(1.37),
            .mass = static_cast<real_t>(183.84),
            .form_factors =
                {
                    static_cast<real_t>(31.5079),
                    static_cast<real_t>(1.629485),
                    static_cast<real_t>(15.682498),
                    static_cast<real_t>(9.446448),
                    static_cast<real_t>(37.960129),
                    static_cast<real_t>(0.000898),
                    static_cast<real_t>(4.885509),
                    static_cast<real_t>(59.980675),
                    static_cast<real_t>(-32.864574),
                },
        },
        {
            .symbol = "Xe",
            .radius = static_cast<real_t>(1.31),
            .mass = static_cast<real_t>(131.29),
            .form_factors =
                {
                    static_cast<real_t>(19.97892),
                    static_cast<real_t>(4.143356),
                    static_cast<real_t>(11.774945),
                    static_cast<real_t>(0.010142),
                    static_cast<real_t>(9.332182),
                    static_cast<real_t>(28.7962),
                    static_cast<real_t>(1.244749),
                    static_cast<real_t>(75.280685),
                    static_cast<real_t>(-6.065902),
                },
        },
        {
            .symbol = "Y",
            .radius = static_cast<real_t>(1.63),
            .mass = static_cast<real_t>(88.906),
            .form_factors =
                {
                    static_cast<real_t>(17.79204),
                    static_cast<real_t>(1.429691),
                    static_cast<real_t>(10.253252),
                    static_cast<real_t>(13.132816),
                    static_cast<real_t>(5.714949),
                    static_cast<real_t>(0.112173),
                    static_cast<real_t>(3.170516),
                    static_cast<real_t>(108.19703),
                    static_cast<real_t>(1.131787),
                },
        },
        {
            .symbol = "Yb",
            .radius = static_cast<real_t>(1.70),
            .mass = static_cast<real_t>(173.05),
            .form_factors =
                {
                    static_cast<real_t>(29.67676),
                    static_cast<real_t>(1.97763),
                    static_cast<real_t>(65.624069),
                    static_cast<real_t>(0.00072),
                    static_cast<real_t>(15.160854),
                    static_cast<real_t>(11.044622),
                    static_cast<real_t>(2.830288),
                    static_cast<real_t>(108.13915),
                    static_cast<real_t>(-60.313812),
                },
        },
        {
            .symbol = "Zn",
            .radius = static_cast<real_t>(1.18),
            .mass = static_cast<real_t>(65.38),
            .form_factors =
                {
                    static_cast<real_t>(14.741002),
                    static_cast<real_t>(3.388232),
                    static_cast<real_t>(6.907748),
                    static_cast<real_t>(0.243315),
                    static_cast<real_t>(4.642337),
                    static_cast<real_t>(11.903689),
                    static_cast<real_t>(2.191766),
                    static_cast<real_t>(63.31213),
                    static_cast<real_t>(-36.915829),
                },
        },
        {
            .symbol = "Zr",
            .radius = static_cast<real_t>(1.54),
            .mass = static_cast<real_t>(91.224),
            .form_factors =
                {
                    static_cast<real_t>(17.859772),
                    static_cast<real_t>(1.310692),
                    static_cast<real_t>(10.911038),
                    static_cast<real_t>(12.319285),
                    static_cast<real_t>(5.821115),
                    static_cast<real_t>(0.104353),
                    static_cast<real_t>(3.512513),
                    static_cast<real_t>(91.777542),
                    static_cast<real_t>(1.124859),
                },
        },
    },
};

constexpr const ElementData *find(std::string_view symbol) noexcept {
  // Binary search implementation for element data retrieval.
  // This function provides O(log N) lookup complexity.
  // @param symbol Chemical symbol.
  // @return Constant pointer to ElementData or nullptr.
  const auto *iterator = std::lower_bound(
      DATABASE.data(), DATABASE.data() + DATABASE.size(), symbol,
      [](const ElementData &data, std::string_view sym) { return data.symbol < sym; });
  if (iterator != DATABASE.data() + DATABASE.size() && iterator->symbol == symbol) {
    return iterator;
  }
  return nullptr;
}
} // namespace detail

} // namespace correlation::physics
