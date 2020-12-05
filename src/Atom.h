#ifndef _ATOM_H_
#define _ATOM_H_

// Minimal structure that represents an atom
struct Atom_Img {
    int                   element_id;
    int                   atom_id;
    std::array<double, 3> position;
};


class Atom {
    // This object represents every atom in the cell
private:
    int number;
    static int NumOfAtoms;
public:
    std::string element;
    int element_id;
    std::array<double, 3> position;
    std::vector<Atom_Img> bonded_atoms;
    // Setters & Getters
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
    // Produce a minimal structure to compute the bond angle.
    Atom_Img GetImage();
    // Get the angle between other atom (Atom_Img) and this object
    double GetAngle(Atom_Img, Atom_Img);
};

#endif /* end of include guard: _ATOM_H_ */
