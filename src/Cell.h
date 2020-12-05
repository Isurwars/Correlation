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
    // Matrix of Bond-lengths
    std::vector<std::vector<double> > bond_length;
    // 3D Tensor of Distances
    std::vector<std::vector<std::vector<double> > > distances;
    // 3D Tensor of Coordination Numbers
    std::vector<std::vector<std::vector<int> > > coordination;
    // 4D Tensor of Angles
    std::vector<std::vector<std::vector<std::vector<double> > > > angles;
    // Matrix of g(r) Histograms
    std::vector<std::vector<double> > g;
    // Matrix of J(r) Histograms
    std::vector<std::vector<double> > J;
    // Matrix of G(r) Histograms
    std::vector<std::vector<double> > G;
    // Matrix of f(theta) Histograms
    std::vector<std::vector<double> > f_theta;
    // Volume of the cell
    double volume;
    // Standar Constructors
    Cell(std::array<double, 6>);
    Cell();
    // Lattice Vectors
    void SetFromVectors(std::vector<double>, std::vector<double>, std::vector<double>);
    void SetLatticeVectors();
    // In Cell corrected positions
    void CorrectPositions();
    void CorrectFracPositions();
    // Populate the Bond length Matrix
    void PopulateBondLength(double);
    //  Read Bond File
    void read_BOND(std::string);
    // RDF calculation (max distance between atoms)
    void RDF(double);
    // Coordination Numbers Calculation
    void Nc();
    // Bond-Angle Calulation ()
    void PAD(bool = true);
    // RDF Histograms  (max distance between atoms, bin width)
    void RDF_Histogram(std::string, double = 20.0, double = 0.05);
    // Nc Histograms ()
    void Nc_Histogram(std::string);
    // PAD Histograms ()
    void PAD_Histogram(std::string, double = 180.0, double = 1.0);
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
