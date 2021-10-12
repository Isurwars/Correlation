#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <list>
#include <array>
#include <vector>
#include <cmath>
#include <map>
#include <regex>
#include <numeric>

#include "Atom.h"
#include "Cell.h"
#include "Constants.h"

/*
 * Generic function to find if an element of any type exists in a vector,
 * if true, then returns the index.
 */
template <typename T>
std::pair<bool, int> findInVector(const std::vector<T> & vecOfElements, const T  & element)
{
    std::pair<bool, int> result;

    // Find given element in vector
    auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

    if (it != vecOfElements.end()) {
        result.second = std::distance(vecOfElements.begin(), it);
        result.first  = true;
    } else {
        result.first  = false;
        result.second = -1;
    }
    return result;
}// findInVector

// Lattice parameters constructor
Cell::Cell(std::array<double, 6> lat)
{
    this->lattice_parameters = lat;
    this->atoms    = { };
    this->elements = { };
    this->SetLatticeVectors();
};


// Default constructor
Cell::Cell()
{
    this->lattice_parameters = { 1.0, 1.0, 1.0, 90, 90, 90 };
    this->atoms    = { };
    this->elements = { };
    this->SetLatticeVectors();
};

// Lattice vector constructor
void Cell::SetFromVectors(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3)
{
    double v1_n, v2_n, v3_n, aux, a, b, c;

    v1_n = sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), 0));
    v2_n = sqrt(std::inner_product(v2.begin(), v2.end(), v2.begin(), 0));
    v3_n = sqrt(std::inner_product(v3.begin(), v3.end(), v3.begin(), 0));
    aux  = std::inner_product(v2.begin(), v2.end(), v3.begin(), 0);
    a    = acos(aux / (v2_n * v3_n)) * constants::rad2deg;
    aux  = std::inner_product(v1.begin(), v1.end(), v3.begin(), 0);
    b    = acos(aux / (v1_n * v3_n)) * constants::rad2deg;
    aux  = std::inner_product(v2.begin(), v2.end(), v1.begin(), 0);
    c    = acos(aux / (v2_n * v1_n)) * constants::rad2deg;
    this->lattice_parameters = { v1_n, v2_n, v3_n, a, b, c };
    this->atoms    = { };
    this->elements = { };
    this->SetLatticeVectors();
};


// Calculate the lattice vectors from the lattice parameters.
void Cell::SetLatticeVectors()
{
    double A     = this->lattice_parameters[0];
    double B     = this->lattice_parameters[1];
    double C     = this->lattice_parameters[2];
    double alpha = this->lattice_parameters[3] * constants::deg2rad;
    double beta  = this->lattice_parameters[4] * constants::deg2rad;
    double gamma = this->lattice_parameters[5] * constants::deg2rad;

    this->v_a_ = { A,
                   0.0,
                   0.0 };
    this->v_b_ = { B * cos(gamma),
                   B * sin(gamma),
                   0.0 };
    this->v_c_ = { C * cos(beta),
                   C * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                   C * sqrt(1
                     - pow(cos(alpha), 2)
                     - pow(cos(beta), 2)
                     - pow(cos(gamma), 2)
                     + 2 * cos(alpha) * cos(beta) * cos(gamma))
                   / sin(gamma) };
    this->volume = (this->v_c_[0] * (this->v_a_[1] * this->v_b_[2] - this->v_a_[2] * this->v_b_[1])
      - this->v_c_[1] * (this->v_a_[0] * this->v_b_[2] - this->v_a_[2] * this->v_b_[0])
      + this->v_c_[2] * (this->v_a_[0] * this->v_b_[1] - this->v_a_[1] * this->v_b_[0]));
}

// Correct the initial positions to in-cell positions
void Cell::CorrectPositions()
{
    double i, j, k;
    int i_, j_, k_, m;
    std::array<double, 3> aux_pos;
    std::list<Atom>::iterator MyAtom;
    // Store the number of atoms per element
    std::vector<int> temp_num_atoms(this->elements.size(), 0);


    for (MyAtom = this->atoms.begin();
      MyAtom != this->atoms.end();
      MyAtom++)
    {
        temp_num_atoms[MyAtom->element_id]++;
        aux_pos = MyAtom->position;
        k       = aux_pos[2] / this->v_c_[2];
        for (m = 0; m < 3; m++) {
            aux_pos[m] -= k * this->v_c_[m];
        }

        j = aux_pos[1] / this->v_b_[1];
        for (m = 0; m < 3; m++) {
            aux_pos[m] -= j * this->v_b_[m];
        }
        i = aux_pos[0] / this->v_a_[0];
        if (i >= -0.000000000000001) {
            i_ = trunc(i);
        } else {
            i_ = trunc(i) - 1;
        }
        if (j >= -0.000000000000001) {
            j_ = trunc(j);
        } else {
            j_ = trunc(j) - 1;
        }
        if (k >= -0.000000000000001) {
            k_ = trunc(k);
        } else {
            k_ = trunc(k) - 1;
        }
        for (m = 0; m < 3; m++) {
            MyAtom->position[m] -= (i_ * this->v_a_[m] + j_ * this->v_b_[m] + k_ * this->v_c_[m]);
        }
    }

    this->element_numbers = temp_num_atoms;
} // Cell::CorrectPositions

// Correct the fractional positions to absolute positions
void Cell::CorrectFracPositions()
{
    double i, j, k;
    std::list<Atom>::iterator MyAtom;

    for (MyAtom = this->atoms.begin();
      MyAtom != this->atoms.end();
      MyAtom++)
    {
        i = MyAtom->position[0];
        j = MyAtom->position[1];
        k = MyAtom->position[2];

        MyAtom->position[0] = (i * this->v_a_[0]
          + j * this->v_b_[0]
          + k * this->v_c_[0]);
        MyAtom->position[1] = (i * this->v_a_[1]
          + j * this->v_b_[1]
          + k * this->v_c_[1]);
        MyAtom->position[2] = (i * this->v_a_[2]
          + j * this->v_b_[2]
          + k * this->v_c_[2]);
    }
} // Cell::CorrectFracPositions

// Populate Bond_length matrix
void Cell::PopulateBondLength(double Bond_Factor)
{
    std::list<Atom>::iterator MyAtom;
    std::pair<bool, int> MyId;
    int i, j;
    double aux;

    // Number of elements in the Cell
    const int n = this->elements.size();
    // Initialize Bond length matrix (nxn) as zeros
    std::vector<std::vector<double> > temp_matrix(n, std::vector<double>(n, 0.0));

    // Iterate in Atoms list to assign the id in the matrix to every atom.
    for (MyAtom = this->atoms.begin();
      MyAtom != this->atoms.end();
      MyAtom++)
    {
        MyId = findInVector(this->elements, MyAtom->element);
        MyAtom->element_id = MyId.second;
    }

    for (i = 0; i < n; i++) {
        aux = Covalent_Radii(this->elements[i]);
        for (j = 0; j < n; j++) {
            temp_matrix[i][j] = (aux + Covalent_Radii(this->elements[j])) * Bond_Factor;
        }
    }
    this->bond_length = temp_matrix;
}// Cell::PopulateBondlength

void Cell::read_BOND(std::string file_name)
{
    /*
     * This function reads the in_bond_file to populate the bond_length Tensor.
     * The file should be in the format:
     * element_B element_A distance(in Angstroms)
     *
     * For example:
     *
     * Si Si 2.29
     * Mg Mg 2.85
     * C  C  1.55
     * C  Si 1.86
     * Si Mg 2.57
     * C  Mg 2.07
     *
     * Any missing pair of elements will use the bond_parameter as a default.
     */
    std::ifstream myfile(file_name);
    std::string line;
    std::smatch match;
    std::pair<bool, int> MyIdA, MyIdB;
    int i, j;
    double dist;

    /*
     * Every line should have two elements and a bond length separeted by spaces:
     *
     * element_A element_B bond_length
     */

    std::regex regex_bond("^([A-Z][a-z]?)"
      "(\\s+)"
      "([A-Z][a-z]?)"
      "(\\s+[-+]?[0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?)");

    if (myfile.is_open()) {
        /* Check if the file is open */
        while (std::getline(myfile, line)) {
            /* Read line by line */
            if (std::regex_search(line, match, regex_bond)) {
                /* Bond found */
                MyIdA = findInVector(this->elements, std::string(match.str(1).data()));
                if (MyIdA.first) i = MyIdA.second;
                MyIdB = findInVector(this->elements, std::string(match.str(3).data()));
                if (MyIdB.first) j = MyIdB.second;
                dist = std::stof(match.str(4).data());
                if (MyIdA.first && MyIdB.first) {
                    this->bond_length[i][j] = dist;
                    this->bond_length[j][i] = dist;
                }
            }
        }
    }
} // read_BOND

// Cell::UpdateProgressBar
void Cell::UpdateProgressBar(double pos)
{
    int barWidth = 50;
    int progress = round(pos * barWidth);
    std::cout << "\r[";

    for (int i = 0; i < barWidth; ++i) {
        if (i < progress) std::cout << "=";
        else if (i == progress) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(pos * 100.0) << " %";
    std::cout.flush();
}// Cell::UpdateProgressBar

// Cell::RDF_MultiThreading
void Cell::RDF_MP(double r_cut)
{
    std::list<Atom>::iterator atom_A;
    std::list<Atom>::iterator atom_B;
    Atom img_atom;
    int id_A, id_B, i, j, k, i_, j_, k_;
    double aux_dist;
    double h_       = 1.0 / this->atoms.size();
    double progress = 0.0;


    // Number of elements in the Cell
    const int n = this->elements.size();
    // This matrix stores the distances between different types of elements,
    // there are at most nxn partials, off-diagonal partials are symmetric.
    std::vector<std::vector<std::vector<double> > > temp_dist(n, std::vector<std::vector<double> >(n,
      std::vector<double>(0)));

    // Correct the atom positions to be inside the cell.
    this->CorrectPositions();

    /*  Our cell is big enought to be divided at least 2 times by the r_cut
     *  in at least one of it's dimensions. This Cell is now our MegaCell,
     *  and we can divide this MegaCell in smaller cells.
     *  We then can calculate the partial distance matrix by only looking
     *  to each individual cell and it's neighbours in a 3x3x3
     *  supercell created by concatenations of the smallers cells.
     */

    k_ = ceil(this->v_c()[2] / r_cut);
    j_ = ceil(this->v_b()[1] / r_cut);
    i_ = ceil(this->v_a()[0] / r_cut);


    /*
     * This is the main loop in RDF, we need to iterate between all atoms
     * in the cell to all the atoms in the supercell created above.
     *
     * This loop is the most time-consuming part of the code, scaling as n^2.
     *
     * We reduce the computation time in half by only iterating for A>B and
     * duplicating the distance because distance from atom_A to atom_B is
     * equal to the distance from atom_B to atom_A.
     *
     * This code can also be improved through paralellization because each
     * iteration is independent of the others, so changing to a for_each is
     * desired to improve this code further.
     */
    for (atom_A = this->atoms.begin();
      atom_A != this->atoms.end();
      atom_A++)
    {
        progress += h_;
        this->UpdateProgressBar(progress);
        id_A = findInVector(this->elements, atom_A->element).second;
        for (atom_B = this->atoms.begin();
          atom_B != this->atoms.end();
          atom_B++)
        {
            // Excluding self-interactions
            if (atom_A->GetNumber() != atom_B->GetNumber()) {
                id_B     = findInVector(this->elements, atom_B->element).second;
                img_atom = *atom_B;
                for (i = -i_; i <= i_; i++) {
                    for (j = -j_; j <= j_; j++) {
                        for (k = -k_; k <= k_; k++) {
                            img_atom.position     = atom_B->position;
                            img_atom.position[0] += i * this->v_a()[0] + j * this->v_b()[0] + k * this->v_c()[0];
                            img_atom.position[1] += j * this->v_b()[1] + k * this->v_c()[1];
                            img_atom.position[2] += k * this->v_c()[2];
                            aux_dist = atom_A->Distance(img_atom);
                            if (aux_dist <= r_cut) {
                                temp_dist[id_A][id_B].push_back(aux_dist);
                                if (aux_dist <= this->bond_length[atom_A->element_id][img_atom.element_id]) {
                                    atom_A->bonded_atoms.push_back(img_atom.GetImage());
                                    if (aux_dist < 0.1) {
                                        std::cout << std::endl << "ERROR: The atoms:" << std::endl
                                                  << atom_A->element << "_" << atom_A->GetNumber()
                                                  << " in position ("
                                                  << atom_A->position[0] << ", "
                                                  << atom_A->position[1] << ", "
                                                  << atom_A->position[2] << ")," << std::endl
                                                  << atom_B->element << "_" << atom_B->GetNumber()
                                                  << " in position ("
                                                  << img_atom.position[0] << ", "
                                                  << img_atom.position[1] << ", "
                                                  << img_atom.position[2] << ")." << std::endl
                                                  << "Have a distance less than 10 pm." << std::endl;
                                        exit(1);
                                    }
                                    if (aux_dist < 0.5) {
                                        std::cout << std::endl << "WARNING: The atoms:" << std::endl
                                                  << atom_A->element << "_" << atom_A->GetNumber()
                                                  << " in position ("
                                                  << atom_A->position[0] << ", "
                                                  << atom_A->position[1] << ", "
                                                  << atom_A->position[2] << ")," << std::endl
                                                  << atom_B->element << "_" << atom_B->GetNumber()
                                                  << " in position ("
                                                  << img_atom.position[0] << ", "
                                                  << img_atom.position[1] << ", "
                                                  << img_atom.position[2] << ")." << std::endl
                                                  << "Have a distance less than the Bohr Radius." << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "\r[==================================================] 100 %" << std::endl;
    this->distances = temp_dist;
}// Cell::RDF_MultiThreading

// RDF Calculation
void Cell::RDF(double r_cut)
{
    std::list<Atom>::iterator atom_A;
    std::list<Atom>::iterator atom_B;
    Atom img_atom;
    int id_A, id_B, i, j, k, i_, j_, k_;
    double aux_dist;
    double h_       = 1.0 / this->atoms.size();
    double progress = 0.0;


    // Number of elements in the Cell
    const int n = this->elements.size();
    // This matrix stores the distances between different types of elements,
    // there are at most nxn partials, off-diagonal partials are symmetric.
    std::vector<std::vector<std::vector<double> > > temp_dist(n, std::vector<std::vector<double> >(n,
      std::vector<double>(0)));

    // Correct the atom positions to be inside the cell.
    this->CorrectPositions();

    //  We need to create a big enough supercell to cover a sphere of radius r_cut.
    k_ = ceil(r_cut / this->v_c()[2]);
    j_ = ceil(r_cut / this->v_b()[1]);
    i_ = ceil(r_cut / this->v_a()[0]);

    /*
     * This is the main loop in RDF, we need to iterate between all atoms
     * in the cell to all the atoms in the supercell created above.
     *
     * This loop is the most time-consuming part of the code, scaling as n^2.
     *
     * We reduce the computation time in half by only iterating for A>B and
     * duplicating the distance because distance from atom_A to atom_B is
     * equal to the distance from atom_B to atom_A.
     *
     * This code can also be improved through paralellization because each
     * iteration is independent of the others, so changing to a for_each is
     * desired to improve this code further.
     */
    for (atom_A = this->atoms.begin();
      atom_A != this->atoms.end();
      atom_A++)
    {
        progress += h_;
        this->UpdateProgressBar(progress);
        id_A = findInVector(this->elements, atom_A->element).second;
        for (atom_B = this->atoms.begin();
          atom_B != this->atoms.end();
          atom_B++)
        {
            // Excluding self-interactions
            if (atom_A->GetNumber() != atom_B->GetNumber()) {
                id_B     = findInVector(this->elements, atom_B->element).second;
                img_atom = *atom_B;
                for (i = -i_; i <= i_; i++) {
                    for (j = -j_; j <= j_; j++) {
                        for (k = -k_; k <= k_; k++) {
                            img_atom.position     = atom_B->position;
                            img_atom.position[0] += i * this->v_a()[0] + j * this->v_b()[0] + k * this->v_c()[0];
                            img_atom.position[1] += j * this->v_b()[1] + k * this->v_c()[1];
                            img_atom.position[2] += k * this->v_c()[2];
                            aux_dist = atom_A->Distance(img_atom);
                            if (aux_dist <= r_cut) {
                                temp_dist[id_A][id_B].push_back(aux_dist);
                                if (aux_dist <= this->bond_length[atom_A->element_id][img_atom.element_id]) {
                                    atom_A->bonded_atoms.push_back(img_atom.GetImage());
                                    if (aux_dist < 0.1) {
                                        std::cout << std::endl << "ERROR: The atoms:" << std::endl
                                                  << atom_A->element << "_" << atom_A->GetNumber()
                                                  << " in position ("
                                                  << atom_A->position[0] << ", "
                                                  << atom_A->position[1] << ", "
                                                  << atom_A->position[2] << ")," << std::endl
                                                  << atom_B->element << "_" << atom_B->GetNumber()
                                                  << " in position ("
                                                  << img_atom.position[0] << ", "
                                                  << img_atom.position[1] << ", "
                                                  << img_atom.position[2] << ")." << std::endl
                                                  << "Have a distance less than 10 pm." << std::endl;
                                        exit(1);
                                    }
                                    if (aux_dist < 0.5) {
                                        std::cout << std::endl << "WARNING: The atoms:" << std::endl
                                                  << atom_A->element << "_" << atom_A->GetNumber()
                                                  << " in position ("
                                                  << atom_A->position[0] << ", "
                                                  << atom_A->position[1] << ", "
                                                  << atom_A->position[2] << ")," << std::endl
                                                  << atom_B->element << "_" << atom_B->GetNumber()
                                                  << " in position ("
                                                  << img_atom.position[0] << ", "
                                                  << img_atom.position[1] << ", "
                                                  << img_atom.position[2] << ")." << std::endl
                                                  << "Have a distance less than the Bohr Radius." << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "\r[==================================================] 100 %" << std::endl;
    this->distances = temp_dist;
} // Cell::RDF

void Cell::CoordinationNumber()
{
    int max_Nc = 0;
    int i;
    std::list<Atom>::iterator MyAtom;
    std::vector<Atom_Img>::iterator atom_A;

    // Search for the maximum number of bonds in the cell
    for (MyAtom = this->atoms.begin();
      MyAtom != this->atoms.end();
      MyAtom++)
    {
        if (max_Nc < int(MyAtom->bonded_atoms.size())) max_Nc = int(MyAtom->bonded_atoms.size());
    }
    const int n = this->elements.size();
    const int m = max_Nc + 2;
    // Create the Tensor nxnxm and initialize it with zeros
    std::vector<std::vector<std::vector<int> > > temp_nc(n, std::vector<std::vector<int> >(n,
      std::vector<int>(m, 0)));
    // Search for the number of bonds per atom per element
    for (MyAtom = this->atoms.begin();
      MyAtom != this->atoms.end();
      MyAtom++)
    {
        std::vector<int> aux(n, 0);
        for (atom_A = MyAtom->bonded_atoms.begin();
          atom_A != MyAtom->bonded_atoms.end();
          atom_A++)
        {
            aux[atom_A->element_id]++;
        }
        for (i = 0; i < n; i++) {
            temp_nc[MyAtom->element_id][i][aux[i]]++;
        }
    }
    this->coordination = temp_nc;
}// Cell::CoordinationNumber

void Cell::PAD(bool degree)
{
    std::list<Atom>::iterator MyAtom;
    std::vector<Atom_Img>::iterator atom_A, atom_B;
    double factor = 1.0;

    if (degree) factor = constants::rad2deg;

    // NxNxN Tensor to store the PAD
    const int n = this->elements.size();
    std::vector<std::vector<std::vector<std::vector<double> > > > temp_pad(n,
      std::vector<std::vector<std::vector<double> > >(n, std::vector<std::vector<double> >(n,
      std::vector<double>(0))));

    /*
     * This is the main loop to calculate the angles formed by every three atoms.
     * The connected atoms are calculated in Cell::RDF, and MUST be called first.
     *
     * The first loop iterates over every atom in the cell.
     * The second and third loop iterate over the connected atoms in every atom instance.
     * These three loops populate a 3D tensor of vectors, the indices of the Tensor
     * represent the three indices of the element in Cell::elements.
     *
     * By default the angle returned is given in degrees.
     */
    for (MyAtom = this->atoms.begin();
      MyAtom != this->atoms.end();
      MyAtom++)
    {
        for (atom_A = MyAtom->bonded_atoms.begin();
          atom_A != MyAtom->bonded_atoms.end();
          atom_A++)
        {
            for (atom_B = MyAtom->bonded_atoms.begin();
              atom_B != MyAtom->bonded_atoms.end();
              atom_B++)
            {
                if (atom_A->atom_id != atom_B->atom_id) {
                    temp_pad[atom_A->element_id][MyAtom->element_id][atom_B->element_id].push_back(MyAtom->GetAngle(*
                      atom_A, *atom_B) * factor);
                }
            }
        }
    }
    this->angles = temp_pad;
}// Cell::PAD

void Cell::RDF_Histogram(std::string filename, double r_cut, double bin_width, bool normalize)
{
    int n_, m_, i, j, col, row;
    int n = this->distances.size();
    std::string header;
    double num_atoms = this->atoms.size();

    m_ = ceil(r_cut / bin_width);
    n_ = n * (n + 1) / 2 + 1;

    /*
     * w_ij is the weighting factor for the partial of G_ij
     * There are to normalization commonly used:
     *     - HHS(Atlas) normalization, commonly used by
     *       experimental scientist.
     *     - Ashcroft - Waseda normalization, commonly
     *       used by theoretical scientist.
     * By default we use Ashcroft normalization:
     * w_ij = c_i * c_j (Product of the numeric concentrations)
     * w_ij = (#atoms_i * #atoms_j)/total_number_atoms^2
     * We offer the option to normalize to HHS by using the
     * -n, --normalize option.
     */
    std::vector<double> temp_w_ij(n_, 1.0);

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            temp_w_ij[i + j + 1] = 2.0 * this->element_numbers[i] * this->element_numbers[j] / ( num_atoms * num_atoms);
            if (i == j) temp_w_ij[i + j + 1] *= 0.5;
        }
    }
    this->w_ij = temp_w_ij;


    std::vector<std::vector<double> > temp_hist(n_, std::vector<double>(m_, 0));
    // Fill the r values of the histogram
    for (i = 0; i < m_; i++) {
        temp_hist[0][i] = (i + 0.5) * bin_width;
    }
    col = 0;
    // Triple loop to iterate over the distances tensor.
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            col++;
            for (std::vector<double>::iterator it = this->distances[i][j].begin();
              it != this->distances[i][j].end(); it++)
            {
                row = floor(*it / bin_width);
                if (row < m_) {
                    temp_hist[col][row]++;
                    if (i != j) temp_hist[col][row]++;
                }
            }
        }
    }

    /*
     * Scale the histograms by the factor: 1 / (#atoms * bin_width)
     */

    double w_factor = num_atoms * bin_width;
    for (i = 1; i < n_; i++) {
        if (normalize) {
            w_factor = num_atoms * bin_width * this->w_ij[i];
        }
        std::transform(temp_hist[i].begin(),
          temp_hist[i].end(),
          temp_hist[i].begin(),
          [&w_factor](auto& c){
            return c / w_factor;
        });
    }

    this->J = temp_hist;

    std::ofstream out_file(filename + "_J.csv");
    std::setprecision(6);
    out_file << std::setw(13) << "r (Å),";
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            header = this->elements[i] + "-" + this->elements[j] + " (1/Å),";
            out_file << std::setw(13) << header;
        }
    }
    out_file << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file << std::setw(11) << temp_hist[j][i] << ",";
        }
        out_file << std::endl;
    }
    out_file.close();

    /*
     * Calculate g(r) with the inverse of the J(r) definition:
     * J(r) = 4 * pi * r^2 * rho_0 * g(r)
     */
    // numeric density
    double rho_0 = num_atoms / this->volume;
    for (col = 1; col < n_; col++) {
        for (row = 1; row < m_; row++) {
            temp_hist[col][row] /= 4 * constants::pi * rho_0 * temp_hist[0][row] * temp_hist[0][row];
        }
    }

    this->g = temp_hist;

    std::ofstream out_file2(filename + "_g.csv");
    std::setprecision(6);
    out_file2 << std::setw(13) << "r (Å),";

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            header = this->elements[i] + "-" + this->elements[j] + " (1/Å),";
            out_file2 << std::setw(13) << header;
            if (normalize) {
                temp_w_ij[i + j + 1] = 1.0;
            }
        }
    }
    out_file2 << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file2 << std::setw(11) << temp_hist[j][i] << ",";
        }
        out_file2 << std::endl;
    }
    out_file2.close();

    /*
     * We calclate G(r) with the definition:
     * G(r) = 4 * pi * r * rho_0 * [g(r) - 1];
     *
     * Weighted partials are calculated with:
     * G_ij = 4 * pi * r * rho_0 * [g_ij(r) - w_ij]
     */

    for (col = 1; col < n_; col++) {
        for (row = 1; row < m_; row++) {
            temp_hist[col][row] = 4 * constants::pi * rho_0 * temp_hist[0][row]
              * (temp_hist[col][row] - temp_w_ij[col]);
        }
    }

    this->G = temp_hist;

    std::ofstream out_file3(filename + "_G_.csv");
    std::setprecision(6);
    out_file3 << std::setw(13) << "r (Å),";
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            header = this->elements[i] + "-" + this->elements[j] + " (1/Å),";
            out_file3 << std::setw(13) << header;
        }
    }
    out_file3 << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file3 << std::setw(11) << temp_hist[j][i] << ",";
        }
        out_file3 << std::endl;
    }
    out_file3.close();
}// Cell::RDF_histogram

void Cell::CoordinationNumber_Histogram(std::string filename)
{
    int n_, m_, i, j, col;
    int n = this->elements.size();
    std::string header;

    n_ = n * n + 1;
    m_ = this->coordination[0][0].size();

    std::vector<std::vector<int> > temp_hist(n_, std::vector<int>(m_, 0));
    // Fill the number of bonds values of the histogram
    for (i = 0; i < m_; i++) {
        temp_hist[0][i] = i;
    }
    col = 1;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            temp_hist[col] = this->coordination[i][j];
            col++;
        }
    }


    std::ofstream out_file2(filename + "_Z.csv");
    std::setprecision(6);
    out_file2 << std::setw(13) << "Number (#),";
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            header = this->elements[i] + " by " + this->elements[j] + ",";
            out_file2 << std::setw(13) << header;
        }
    }
    out_file2 << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file2 << std::setw(12) << temp_hist[j][i] << ",";
        }
        out_file2 << std::endl;
    }
    out_file2.close();
}// Cell::CoordinationNumber_histogram

void Cell::SQ(std::string filename, double q_bin_width, double bin_width, double r_cut, bool normalize)
{
    int n_, m_, i, j, row, col;
    int n = this->elements.size();
    std::string header;
    double Trapz     = 0.0;
    double num_atoms = this->atoms.size();

    n_ = this->G.size();
    m_ = this->G[0].size();
    std::vector<double> temp_w_ij(n_, 1.0);
    if (!normalize) {
        temp_w_ij = this->w_ij;
    }

    /*
     * The structure factor S(q) is calculated with:
     * S(q) = 1 + 4*pi*rho_0*(q^-1)*\int{ dr r*sin(qr)*[g(r) - 1]
     * or S(q) = 1 + (q^{-1})*\int{dr sin(q * r) * G(r)}
     */
    std::vector<std::vector<double> > temp_S(n_, std::vector<double>(m_, 0));

    // Fill the q values of the histogram
    for (i = 0; i < m_; i++) {
        temp_S[0][i] = 0.0 + (i + 0.0) * q_bin_width;
    }
    double rho_0 = this->atoms.size() / this->volume;
    double aux   = 4 * constants::pi * rho_0;

    /*
     * Double loop on q (rows) and G_ij (cols)
     */
    for (col = 1; col < n_; col++) {
        for (row = 0; row < m_; row++) {
            /*
             * Integration with Trapezoidal_rule
             */
            Trapz = 0.0;
            if (temp_S[0][row] < 0.2) {
                for (i = 1; i < m_; i++) {
                    Trapz += (1 - 0.5 * pow(temp_S[0][row] * this->g[0][i], 2)) * this->g[col][i];
                }
                temp_S[col][row] = temp_w_ij[col] + Trapz * bin_width * rho_0;
            } else {
                for (i = 1; i < (m_ - 1); i++) {
                    Trapz += std::sin(temp_S[0][row] * this->G[0][i]) * this->G[col][i];
                }
                Trapz += 0.5 * std::sin(temp_S[0][row] * this->G[0][m_ - 1]) * this->G[col][m_ - 1];
                temp_S[col][row] = temp_w_ij[col] + Trapz * (bin_width / temp_S[0][row]);
            }
        }
    }

    this->S = temp_S;

    std::ofstream out_file(filename + "_S_q.csv");
    std::setprecision(6);
    out_file << std::setw(13) << "q (1/Å),";
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            header = this->elements[i] + "-" + this->elements[j] + " (Å),";
            out_file << std::setw(13) << header;
        }
    }
    out_file << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file << std::setw(11) << temp_S[j][i] << ",";
        }
        out_file << std::endl;
    }
    out_file.close();
}// Cell::SQ

void Cell::XRD(std::string filename, double lambda, double theta_min, double theta_max, double bin_width)
{
    int n_, m_, i, j, k, col, row;
    double norm, aux;
    int n = this->elements.size();
    std::string header;

    m_ = ceil((theta_max - theta_min) / bin_width);
    n_ = n * (n + 1) / 2 + 1;

    std::vector<std::vector<double> > temp_hist(n_, std::vector<double>(m_, 0));
    // Fill the r values of the histogram
    for (i = 0; i < m_; i++) {
        temp_hist[0][i] = theta_min + i * bin_width;
    }
    col = 0;
    // Triple loop to iterate over the distances tensor.
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            col++;
            for (std::vector<double>::iterator it = this->distances[i][j].begin();
              it != this->distances[i][j].end(); it++)
            {
                /*
                 * Bragg's Law
                 * n * lambda = 2d sin(theta)
                 * theta = asin((n * lambda) / (2 * d))
                 */
                aux = lambda / (*it * 2.0);
                k   = 1;
                std::cout << "d:" << *it << std::endl;
                while ((k * aux) <= 1.0) {
                    row = floor(((2 * constants::rad2deg * asin(k * aux)) - theta_min) / bin_width);
                    k++;
                    if (0 <= row && row < m_) {
                        temp_hist[col][row]++;
                        if (i != j) temp_hist[col][row]++;
                    }
                }
            }
        }
    }

    // Double loop to find normalization factor
    norm = 0.0;
    for (i = 1; i < n_; i++) {
        for (j = 0; j < m_; j++) {
            norm = std::max(temp_hist[i][j], norm);
        }
    }
    norm *= 0.01;
    // Double loop to normalize PAD_Histogram
    for (i = 1; i < n_; i++) {
        for (j = 0; j < m_; j++) {
            temp_hist[i][j] /= norm;
        }
    }

    this->X = temp_hist;
    std::ofstream out_file(filename + "_XRD.csv");
    std::setprecision(5);
    out_file << std::setw(13) << "2-theta (°),";
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            header = this->elements[i] + "-" + this->elements[j] + " (%),";
            out_file << std::setw(13) << header;
        }
    }
    out_file << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file << std::setw(12) << temp_hist[j][i] << ",";
        }
        out_file << std::endl;
    }

    out_file.close();
}// Cell:XRD

void Cell::PAD_Histogram(std::string filename, double theta_cut, double bin_width)
{
    int n_, m_, i, j, k, h, col, row;
    double norm;
    int n = this->elements.size();
    std::string header;

    /*
     * The number of columns in the output file is:
     *         n       x      (n+1)! / [2 x (n-1)!]
     * Central atom        Combination with repetition
     *                         of n in groups of 2
     * it's reduced to n x (n+1) x n /2
     * one extra column is added for Theta (angle)
     */

    n_ = 1 + (n * n * (n + 1) / 2);
    // from 0 to theta_cut degrees rows
    m_ = std::round(theta_cut / bin_width);
    // Matrix to store the Histograms n_ columns, m_ rows
    std::vector<std::vector<double> > temp_hist(n_, std::vector<double>(m_, 0));
    // Fill the theta values of the histogram
    for (i = 0; i < m_; i++) {
        temp_hist[0][i] = (i + 0.5) * bin_width;
    }
    col = 0;
    // Quadruple loop to iterate over the 3D angle tensor.
    for (i = 0; i < n; i++) {         // i iterates over all central atoms
        for (j = 0; j < n; j++) {     // j iterates over all initial atoms
            for (k = j; k < n; k++) { // k iterates only over half + 1 of the spectrum
                col++;
                for (std::vector<double>::iterator it = this->angles[j][i][k].begin();
                  it != this->angles[j][i][k].end(); it++)
                {
                    row = floor(*it / bin_width);
                    if (row < m_) {
                        temp_hist[col][row]++;
                    }
                }
                // Remove double count when j == k
                if (j == k) {
                    for (h = 0; h < m_; h++) {
                        temp_hist[col][h] /= 2.0;
                    }
                }
            }
        }
    }
    // Double loop to find normalization factor
    norm = 0.0;
    for (i = 1; i < n_; i++) {
        for (j = 0; j < m_; j++) {
            norm += temp_hist[i][j];
        }
    }
    norm *= bin_width;

    // Double loop to normalize PAD_Histogram
    for (i = 1; i < n_; i++) {
        for (j = 0; j < m_; j++) {
            temp_hist[i][j] /= norm;
        }
    }

    this->f_theta = temp_hist;
    std::ofstream out_file(filename + "_PAD.csv");
    std::setprecision(5);
    out_file << std::setw(13) << "theta (°),";
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = j; k < n; k++) {
                header = this->elements[j] + "-" + this->elements[i] + "-" + this->elements[k] + ",";
                out_file << std::setw(12) << header;
            }
        }
    }
    out_file << std::endl << std::fixed;

    for (i = 0; i < m_; i++) {
        for (j = 0; j < n_; j++) {
            out_file << std::setw(11) << temp_hist[j][i] << ",";
        }
        out_file << std::endl;
    }

    out_file.close();
}// Cell::PAD_Histogram
