#ifndef SRC_ATOM_H_
#define SRC_ATOM_H_

/* ---------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2021 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * ----------------------------------------------------------------------
 */
#include <array>
#include <string>
#include <vector>

// Minimal structure that represents an atom
struct Atom_Img {
  int                   element_id;
  int                   atom_id;
  std::array<double, 3> position;
};


class Atom {
  // This object represents every atom in the cell
 private:
  int id;
  static int NumOfAtoms;

 public:
  std::string element;
  int element_id;
  std::array<double, 3> position;
  std::vector<Atom_Img> bonded_atoms;
  std::vector<Atom_Img> second_shell_atoms;

  // Setters & Getters
  int GetID() { return id; }
  void SetID(int num) { this->id = num; }

  // Constructors
  Atom(std::string, std::array<double, 3>);
  Atom();
  void SetAll(std::string, std::array<double, 3>);

  // Default functions for the object Atom
  static int GetNumberOfAtoms() { return NumOfAtoms; }
  double Distance(const Atom&);

  // Produce a minimal structure to compute the bond angle.
  Atom_Img GetImage();

  // Get the angle between other atom (Atom_Img) and this object
  double GetAngle(Atom_Img, Atom_Img);

  // Get bonded atoms id
  std::vector<int> GetBondedAtomsID();
};
#endif  // SRC_ATOM_H_
