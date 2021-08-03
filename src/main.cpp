// This program calculate correlation functions for a group of atoms
#include <iostream> // Standar IO library
#include <list>     // To handle the list of atoms
#include <array>    // Handle the array of parameters
#include <vector>   // Most of the output is stored in vectors
#include <getopt.h> // For the argument parsing
#include <string>   // String manipulation and algorithms
#include <algorithm>// for_each

#include "Atom.h"      // Atom Class
#include "Cell.h"      // Cell Class
#include "ReadFiles.h" // File reading and parsing

/*
 * Currenly filesystem is not fully implemented in gcc 10.1.
 * This is an ugly hack to try to make this code compatible
 * in Linux, Windows and Mac.
 * THIS SHOULD HAVE BEEN IMPLEMENTED SINCE 2017!!!, NOT A
 * SINGLE COMPILER HAS ADDED FILESYSTEM SUPPORT YET!!!.
 * Hope this is no longer requiered in the future.
 * #include <filesystem>
 * namespace fs = std::filesystem;
 */
std::pair<std::string, std::string> GetExtension(std::string filename)
{
    std::pair<std::string, std::string> result;
    size_t i = filename.rfind('.', filename.length());

    if (i != std::string::npos) {
        result.first  = filename.substr(0, i);
        result.second = filename.substr(i, filename.length() - 1);
        return result;
    }

    result.first  = filename;
    result.second = "";
    return result;
}

std::string _in_file_name_   = "";
std::string _out_file_name_  = "";
std::string _bond_file_name_ = "";
bool _bond_in_file_  = false;
double _r_cut_       = 20.0;
double _bin_w_       = 0.05;
double _q_bin_w_     = 0.15707963;
double _bond_par_    = 1.3;
double _angle_bin_w_ = 1.0;

void PrintHelp()
{
    const std::string help_txt = "CORRELATION\n"
      "DESCRIPTION:\n"
      "  This program calculates the main correlation functions of a material:\n"
      "    - Radial Distribution Function (J(r)).\n"
      "    - Pair Distribution Function (g(r)).\n"
      "    - Coordination Number (CN).\n"
      "    - Plane-Angle Distribution (PAD).\n\n"
      "USAGE:\n"
      "  The minimal argument is a structure file, this program requires a file\n"
      "  that contains atom positions, crystal structure and composition.\n"
      "  Supported structure files are:\n"
      "    -*.CAR   Materials Studio structure file.\n"
      "    -*.CELL  CASTEP structure file.\n"
      "    -*.dat   ONETEP structure file.\n\n"
      "OPTIONS:\n"
      "    HELP OPTIONS:"
      "      -h, --help\n"
      "        Display this help text.\n\n"
      "    RADIAL OPTIONS:\n"
      "      -r, --r_cut\n"
      "        Cutoff radius in the calculation of g(r), G(r) and J(r). The default\n"
      "        radius it's set to 2 nm. The maximum recomended radius is the same as\n"
      "        shortest length of the periodic boundary conditions (PBC), anything\n"
      "        above this PBC value can be affected by periodic interactions.\n\n"
      "      -w, --bin_width\n"
      "        Width of the histograms for g(r) and J(r), the default is 0.05 nm.\n\n"
      "    BOND-ANGLE OPTIONS:\n"
      "      -a, --angle_bin_width\n"
      "        Width of the histograms for the PAD, default set to 1.0°.\n\n"
      "      -b, --bond_parameter\n"
      "        The ideal covalent bond length is the sum of covalent radii\n"
      "        of the two atoms. The criterion used to consider atoms as bonded\n"
      "        is the following:\n"
      "            0.6 * Sum_radii < distance < bond_parameter * Sum_radii.\n"
      "        By default the bond_parameter is set to 1.30, as a rule of thumb.\n"
      "        The default should work for most crystalline materials,\n"
      "        as well as most covalent non-crystalline materials.\n"
      "        For amorphous and liquid materials the bond_parameter should be\n"
      "        increased to match the desired distance to cut_off the bonds. \n"
      "        A bond_parameter of 1.42 is recomended for amorphous materials.\n\n"
      "      -i, --in_bond_file"
      "        The input file with the bond distances for every pair of elements\n"
      "        in the corresponding input structure. The file should have the\n"
      "        following format:\n"
      "             Si Si 2.29\n"
      "             Mg Mg 2.85\n"
      "             C  C  1.55\n"
      "             C  Si 1.86\n"
      "             Si Mg 2.57\n"
      "             C  Mg 2.07\n"
      "        If any of the pairs is missing in the input file, the corresponding\n"
      "        bond distance will be set using the bond_parameter(1.30 by default).\n\n"
      "    OUTPUT OPTIONS:\n"
      "      -o, --out_file\n"
      "        The output file name, by default the input seed name will be used.\n\n"
      "CREATOR:\n"
      "  This program was created by PhD. Isaias Rodriguez Aguirre, November 2020.\n"
      "  e-mail: isurwars@gmail.com\n"
      "ACKNOWLEDGMENTS:\n"
      "  This software was created during a Posdoctoral fellowship in IIM-UNAM.\n"
      "  Thanks to DGAPA-UNAM for the financial support during my fellowship.\n"
      "  And their continious support as part of the Workgroup projects with IDs:\n"
      "  IN104617 and IN116520.\n";
    std::cout << help_txt;

    exit(1);
} // PrintHelp

void ArgParser(int argc, char ** argv)
{
    const char * const short_opts = "a:b:hi:o:r:w:";
    const option long_opts[]      = {
        { "angle_bin_width", required_argument, nullptr, 'a' },
        { "bond_parameter",  required_argument, nullptr, 'b' },
        { "help",            no_argument,       nullptr, 'h' },
        { "in_bond_file",    required_argument, nullptr, 'i' },
        { "out_file",        required_argument, nullptr, 'o' },
        { "r_cut",           required_argument, nullptr, 'r' },
        { "bin_width",       required_argument, nullptr, 'w' },
        { nullptr,           no_argument,       nullptr, 0   }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (opt == -1)
            break;

        switch (opt) {
            case 'a': // -a or --angle_bin_width
                try{
                    _angle_bin_w_ = std::stof(optarg);
                }
                catch (const std::exception& e) {
                    std::cout << "Invalid input argument: '"
                              << optarg
                              << "' in angle_bin_width parameter '-a' "
                              << "(Real number expected)."
                              << std::endl;
                    exit(1);
                }
                break;
            case 'b': // -b or --bond_parameter
                try{
                    _bond_par_ = std::stof(optarg);
                }
                catch (const std::exception& e) {
                    std::cout << "Invalid input argument: '"
                              << optarg
                              << "' in bond_parameter '-b' "
                              << "(Real number expected)."
                              << std::endl;
                    exit(1);
                }
                break;
            case 'h': // -h or --help
                PrintHelp();
                break;
            case 'i': // -i or --in_bond_file
                _bond_in_file_   = true;
                _bond_file_name_ = optarg;
                break;
            case 'o': // -o or --out_file
                _out_file_name_ = optarg;
                break;
            case 'r': // -r ot --r_cut
                try{
                    _r_cut_ = std::stof(optarg);
                }
                catch (const std::exception& e) {
                    std::cout << "Invalid input argument: '"
                              << optarg
                              << "' in r_cut parameter '-r' "
                              << "(Real number expected)."
                              << std::endl;
                    exit(1);
                }
                break;
            case 'w': // -w or -- bin_width
                try{
                    _bin_w_ = std::stof(optarg);
                }
                catch (const std::exception& e) {
                    std::cout << "Invalid input argument: '"
                              << optarg
                              << "' in bin_width parameter '-w' "
                              << "(Real number expected)."
                              << std::endl;
                    exit(1);
                }
                break;
            case '?': // Unrecognized option
                /* getopt_long already printed an error message. */
                exit(1);
                break;

            default:
                break;
        }
    }
    if (optind < argc - 1) {
        std::cout << "Only one structure file most be provided." << std::endl;
        exit(1);
    } else {
        if (argc > 1) {
            _in_file_name_ = argv[optind];
        } else {
            PrintHelp();
        }
    }
} // ArgParser

int main(int argc, char ** argv)
{
    Cell MyCell; // Struture to be analized
    std::list<Atom>::iterator MyAtom;
    std::pair<std::string, std::string> file_ext;

    /*
     * Parse the input argumets provided
     */
    ArgParser(argc, argv);

    /*
     * Read the structure file
     */
    file_ext = GetExtension(_in_file_name_);
    if (_out_file_name_ == "") _out_file_name_ = file_ext.first;
    std::string MyExt = file_ext.second;
    // Convert extension back to lower case
    std::for_each(MyExt.begin(), MyExt.end(), [](char & c){
        c = ::tolower(c);
    });
    if (MyExt == ".car") {
        MyCell = read_CAR(_in_file_name_);
    } else if (MyExt == ".cell") {
        MyCell = read_CELL(_in_file_name_);
    } else if (MyExt == ".dat") {
        MyCell = read_ONETEP_DAT(_in_file_name_);
    } else {
        std::cout << "File: " << MyExt << " currently not supported." << std::endl;
        PrintHelp();
    }
    std::cout << "File " << _in_file_name_ << " opened successfully." << std::endl;

    // Create Bond distance Matrix and element_ids
    MyCell.PopulateBondLength(_bond_par_);
    if (_bond_in_file_) {
        /*Read the Bond Distances from the external file*/
        MyCell.read_BOND(_bond_file_name_);
    }

    /*
     * This function calculates the distances between every pair of atoms
     * in the structure. A supercell method is used to create the images
     * in all directions up to _r_cut_ distance.
     */
    MyCell.RDF(_r_cut_);

    /*
     * This function calculates the partial coordination number for pairs of
     * elements. Bonded Atoms use the same parameters for PAD.
     */
    MyCell.Nc();

    /*
     * This function calculates the angle between every atom and all pairs
     * of bonded atoms. The bonded atoms are calculated in Cell::RDF and it
     * must be called first.
     */
    MyCell.PAD();

    /*
     * This function uses the distances to calculate g(r), G(r) and J(r).
     * The _r_cut_ parameter is the cutoff distance,
     * the _bin_w_ parameter is the bin width to be used.
     */
    std::cout << "Writing output files: " << _out_file_name_ << std::endl;
    MyCell.RDF_Histogram(_out_file_name_, _r_cut_, _bin_w_);
    MyCell.Nc_Histogram(_out_file_name_);

    /*
     * This function calculate S(Q) as the Fourier Transform of G(r).
     * The _r_cut_ parameter is the cutoff distance in r-space,
     * the _q_bin_w_ parameter is the bin width to be used in q-space.
     */
    MyCell.SQ(_out_file_name_, _r_cut_, _q_bin_w_);

    /*
     * This function uses the angles to calculate the PAD.
     * The theta_max parameter is the maximum angle to compute,
     * the _bin_w_ parameter is the bin width to be used.
     */
    MyCell.PAD_Histogram(_out_file_name_, 180.0, _angle_bin_w_);

    std::cout << "Job in " << _in_file_name_ << " finished successfully." << '\n';
    return 0;
} // main
