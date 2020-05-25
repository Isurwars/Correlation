#ifndef _CELL_H_
#define _CELL_H_

class Cell {
    // This object contains the lattice parameters,
    // and the atoms that compose the material
private:
    // Cell Vectors private modification, public read with getter
    std::array<double, 3> v_a_;
    std::array<double, 3> v_b_;
    std::array<double, 3> v_c_;
public:
    // Lattice Parameters
    std::array<double, 6> lattice_parameters;
    // List of atoms in the Cell
    std::list<Atom> atoms;
    // List of elements in Cell
    std::vector<std::string> elements;
    // Tensor of Distances
    std::vector<std::vector<std::vector<double> > > Distances;
    // Matrix of J(r) Histograms
    std::vector<std::vector<double> > J;
    // Matrix of g(r) Histograms
    std::vector<std::vector<double> > g;
    // Volume of the cell
    double volume;
    // Standar Constructors
    Cell(std::array<double, 6>);
    Cell();
    // Lattice Vectors
    void SetFromVectors(std::array<double, 3>, std::array<double, 3>, std::array<double, 3>);
    void SetLatticeVectors();
    // In Cell corrected positions
    void CorrectPositions();
    void CorrectFracPositions();
    // RDF calculation (max distamce between atoms)
    void RDF(double);
    // RDF Histograms  (max distance between atoms, bin width)
    void Histogram(double, double, std::string);
    // Read Only Lattive Vectors
    std::array<double, 3> v_a()
    {
        return v_a_;
    };
    std::array<double, 3> v_b()
    {
        return v_b_;
    };
    std::array<double, 3> v_c()
    {
        return v_c_;
    };
};

#endif /* end of include guard: _CELL_H_ */
