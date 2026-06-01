// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "app/AppBackend.hpp"
#include <string>
#include <vector>
#include <filesystem>

namespace correlation::app {

/**
 * @struct Preset
 * @brief Represents a named set of user options for materials simulation analysis.
 */
struct Preset {
    std::string name;
    std::string description;
    ProgramOptions options;
};

/**
 * @class PresetManager
 * @brief Manages saving and loading analysis parameter presets to/from JSON files.
 */
class PresetManager {
public:
    /**
     * @brief Gets the folder path where user presets are saved (~/.correlation/presets/).
     */
    static std::filesystem::path presetsDirectory();

    /**
     * @brief Loads all saved presets from the presets directory.
     */
    static std::vector<Preset> loadAll();

    /**
     * @brief Saves a preset to a JSON file in the presets directory.
     */
    static void save(const Preset &preset);

    /**
     * @brief Deletes a saved preset by name.
     */
    static void remove(const std::string &name);

    /**
     * @brief Serializes a Preset object to a JSON string.
     */
    static std::string toJson(const Preset &preset);

    /**
     * @brief Deserializes a Preset object from a JSON string.
     */
    static Preset fromJson(const std::string &json);
};

} // namespace correlation::app
