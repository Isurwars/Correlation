#include <iostream>
#include <list>
#include <array>
#include <vector>
#include <cmath>
#include <string>

#include "Atom.h"
#include "Cell.h"
// All of functions of Class Atom
int Atom::NumOfAtoms = 0;

void Atom::SetAll(std::string ele, std::array<double, 3> pos)
{
    this->element  = ele;
    this->position = pos;
}

// Complete constructor for the object Atom
Atom::Atom(std::string ele, std::array<double, 3> pos)
{
    this->element  = ele;
    this->position = pos;
    this->number   = Atom::NumOfAtoms;
    Atom::NumOfAtoms++;
}

// Default constructor for the object Atom
Atom::Atom()
{
    this->element  = "H";
    this->position = { 0.0, 0.0, 0.0 };
    this->number   = Atom::NumOfAtoms;
    Atom::NumOfAtoms++;
}

// Distance to other atom
double Atom::Distance(const Atom& other_atom)
{
    return sqrt(pow(this->position[0] - other_atom.position[0], 2)
             + pow(this->position[1] - other_atom.position[1], 2)
             + pow(this->position[2] - other_atom.position[2], 2));
}
