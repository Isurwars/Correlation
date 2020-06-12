#include <iostream>
#include <list>
#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>

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

// GetImage
Atom_Img Atom::GetImage()
{
    Atom_Img temp_img;

    temp_img.element_id = this->element_id;
    temp_img.position   = this->position;
    return temp_img;
}

// Get Bond Angle
double Atom::GetAngle(Atom_Img atom_A, Atom_Img atom_B)
{
    std::vector<double> vA, vB;
    double vA_, vB_, aux;

    vA = { atom_A.position[0] - this->position[0],
           atom_A.position[1] - this->position[1],
           atom_A.position[2] - this->position[2] };
    vB = { atom_B.position[0] - this->position[0],
           atom_B.position[1] - this->position[1],
           atom_B.position[2] - this->position[2] };
    vA_  = sqrt(std::inner_product(vA.begin(), vA.end(), vA.begin(), 0));
    vB_  = sqrt(std::inner_product(vB.begin(), vB.end(), vB.begin(), 0));
    aux  = std::inner_product(vA.begin(), vA.end(), vB.begin(), 0);
    aux /= vA_ * vB_;
    if (aux > 1.0) aux = 1.0;
    if (aux < -1.0) aux = -1.0;
    return acos(aux);
}
