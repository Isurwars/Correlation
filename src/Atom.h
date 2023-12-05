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
  /* ---------------------------------------------------------------------
  * This object represents every atom in the cell.
  *
  * The atributes consist of:
  * ID (Unique Identifier),
  * Element (string and id),
  * Postion (array of three doubles),
  * Bonded Atoms (Images of the bonded atoms)
  *
  * As well as several methods to calculate distance between pairs of atoms,
  * and angle between terns of atoms
  * ---------------------------------------------
  */
 private:
  int id;
  std::array<double, 3> _position_;
  std::vector<Atom_Img> _bonded_atoms_;
  int _element_id_;
  std::string _element_;
  static int NumOfAtoms;

 public:
  std::vector<Atom_Img> second_shell_atoms;

  //-------------------------------------------------------------//
  //---------------------- Constructors  ------------------------//
  //-------------------------------------------------------------//

  Atom(std::string, std::array<double, 3>);
  Atom();
  void SetAll(std::string, std::array<double, 3>);

  //-------------------------------------------------------------//
  //-------------------- Setters & Getters ----------------------//
  //-------------------------------------------------------------//

  int GetID() { return this->id; }
  void SetID(int num) { this->id = num; }

  std::array<double, 3> position() { return this->_position_; }
  void SetPosition(std::array<double, 3> pos) { this->_position_ = pos; }

  std::vector<Atom_Img> bonded_atoms() { return this->_bonded_atoms_; }
  std::vector<int> GetBondedAtomsID();
  void AddBondedAtom(Atom_Img);

  int element_id() { return this->_element_id_; }
  void SetElementId(int ele_id) { this->_element_id_ = ele_id; }

  std::string element() { return this->_element_; }
  void SetElement(std::string ele) { this->_element_ = ele; }

  static int GetNumberOfAtoms() { return NumOfAtoms; }


  // Default functions for the object Atom
  double Distance(Atom&);

  // Produce a minimal structure to compute the bond angle.
  Atom_Img GetImage();

  // Get the angle between other atom (Atom_Img) and this object
  double GetAngle(Atom_Img, Atom_Img);
};
#endif  // SRC_ATOM_H_
