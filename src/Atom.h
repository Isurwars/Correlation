#ifndef _ATOM_H_
#define _ATOM_H_


class Atom {
    // This object represent every atom in the group of atoms
private:
    int number;
    std::list<double> bond_distances;
    std::list<double> bonded_atoms;

    static int NumOfAtoms;
public:
    std::string element;
    std::array<double, 3> position;
    // Seters & Geters
    int GetNumber()
    {
        return number;
    };
    void SetNumber(int num)
    {
        this->number = num;
    };
    // Constructors
    Atom(std::string, std::array<double, 3>);
    Atom();
    void SetAll(std::string, std::array<double, 3>);
    // Default functions for the object Atom
    static int GetNumberOfAtoms()
    {
        return NumOfAtoms;
    };
    double Distance(const Atom&);
};

#endif /* end of include guard: _ATOM_H_ */
