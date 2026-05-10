/**
 * @file PdbReader.cpp
 * @brief Implementation of the PDB (.pdb) file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "readers/PdbReader.hpp"
#include "readers/ReaderFactory.hpp"
#include "math/Vector3.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {
bool registered = ReaderFactory::instance().registerReader(
    std::make_unique<PdbReader>());
}

correlation::core::Cell PdbReader::readStructure(
    const std::string& filename,
    std::function<void(float, const std::string&)> progress_callback) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    correlation::core::Cell cell;
    std::string line;
    bool has_box = false;

    while (std::getline(file, line)) {
        if (line.substr(0, 6) == "CRYST1") {
            // CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1
            try {
                double a = std::stod(line.substr(6, 9));
                double b = std::stod(line.substr(15, 9));
                double c = std::stod(line.substr(24, 9));
                cell.setBox(a, b, c);
                has_box = true;
            } catch (...) {}
        } else if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
            // ATOM      1  N   ALA A   1      11.104   6.134  -1.267  1.00  0.00           N
            try {
                std::string symbol = "";
                if (line.length() >= 78) {
                    symbol = line.substr(76, 2);
                    symbol.erase(0, symbol.find_first_not_of(" "));
                    symbol.erase(symbol.find_last_not_of(" ") + 1);
                }
                
                if (symbol.empty()) {
                    // Fallback to atom name columns 12-16
                    std::string atom_name = line.substr(12, 4);
                    atom_name.erase(0, atom_name.find_first_not_of(" "));
                    atom_name.erase(atom_name.find_last_not_of(" ") + 1);
                    
                    for (char c : atom_name) {
                        if (std::isalpha(c)) symbol += c;
                        else break;
                    }
                    // Common PDB fixes (e.g. 1H -> H)
                    if (symbol.length() > 1 && std::isdigit(atom_name[0])) {
                         symbol = atom_name.substr(1, 1);
                    }
                }

                double x = std::stod(line.substr(30, 8));
                double y = std::stod(line.substr(38, 8));
                double z = std::stod(line.substr(46, 8));

                cell.addAtom(correlation::core::Atom(symbol, correlation::math::Vector3<double>(x, y, z)));
            } catch (...) {}
        } else if (line.substr(0, 6) == "ENDMDL" || line.substr(0, 3) == "END") {
            break; 
        }
    }

    if (progress_callback) progress_callback(1.0f, "PDB structure loaded.");
    return cell;
}

correlation::core::Trajectory PdbReader::readTrajectory(
    const std::string& filename,
    std::function<void(float, const std::string&)> progress_callback) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<correlation::core::Cell> frames;
    std::string line;
    correlation::core::Cell current_cell;
    bool in_model = false;
    bool has_box = false;
    double a=0, b=0, c=0;

    while (std::getline(file, line)) {
        if (line.substr(0, 6) == "CRYST1") {
            try {
                a = std::stod(line.substr(6, 9));
                b = std::stod(line.substr(15, 9));
                c = std::stod(line.substr(24, 9));
                has_box = true;
            } catch (...) {}
        } else if (line.substr(0, 5) == "MODEL") {
            current_cell = correlation::core::Cell();
            if (has_box) current_cell.setBox(a, b, c);
            in_model = true;
        } else if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
             try {
                std::string symbol = "";
                if (line.length() >= 78) {
                    symbol = line.substr(76, 2);
                    symbol.erase(0, symbol.find_first_not_of(" "));
                    symbol.erase(symbol.find_last_not_of(" ") + 1);
                }
                if (symbol.empty()) {
                    std::string atom_name = line.substr(12, 4);
                    atom_name.erase(0, atom_name.find_first_not_of(" "));
                    atom_name.erase(atom_name.find_last_not_of(" ") + 1);
                    for (char c : atom_name) {
                        if (std::isalpha(c)) symbol += c;
                        else break;
                    }
                }
                double x = std::stod(line.substr(30, 8));
                double y = std::stod(line.substr(38, 8));
                double z = std::stod(line.substr(46, 8));
                current_cell.addAtom(correlation::core::Atom(symbol, correlation::math::Vector3<double>(x, y, z)));
            } catch (...) {}
        } else if (line.substr(0, 6) == "ENDMDL") {
            frames.push_back(std::move(current_cell));
            in_model = false;
        }
    }

    // If no MODEL tags were found, but we have atoms, it's a single frame PDB
    if (frames.empty()) {
        file.clear();
        file.seekg(0);
        frames.push_back(readStructure(filename, progress_callback));
    }

    if (progress_callback) progress_callback(1.0f, "PDB trajectory loaded.");
    return correlation::core::Trajectory(frames, 1.0);
}

} // namespace correlation::readers
