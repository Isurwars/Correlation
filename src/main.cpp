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
 * Hope this is not requiered in the future.
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

std::string _in_file_name_  = "";
std::string _out_file_name_ = "";
double _r_cut_    = 20.0;
double _bin_w_    = 0.05;
double _bond_len_ = 1.42;

void PrintHelp()
{
    const std::string help_txt = "CORRELATION\n"
      "DESCRIPTION:\n"
      "This program calculates the main correlation functions of a material:\n"
      "- Radial Distribution Function (J(r)).\n"
      "- Pair Distribution Function (g(r)).\n"
      "- Bond Angle Distribution (BAD).\n\n"
      "USAGE:\n"
      "The minimal argument is a structure file, this program requires a file\n"
      "that contains atom positions, crystal structure and composition.\n"
      " Valid structure files are:\n"
      "- *.CAR   Materials Studio structure file.\n"
      "- *.CELL  CASTEP structure file.\n"
      "- *.CIF   Crystallographic Information File.\n\n"
      "OPTIONS:\n"
      "-h, --help\n"
      "Display this help text.\n\n"
      "-r, --r_cut\n"
      "Maximum radios to calculate J(r) and g(r), by default it's set to 2 nm.\n"
      "The maximum recomended radius is the smallest of the periodic boundary\n"
      "conditions (PBC), anything above the smallest PBC can be affected by\n"
      "periodic interactions.\n\n"
      "-w, --bin_width\n"
      "Width of the histograms bins, default set to 0.005 nm.\n\n"
      "-b, --bond_length\n"
      "The ideal bond length is the sum of covalent radii of the two atoms.\n"
      "Our criteria is the following:\n"
      "0.6 * Sum_radii < distance < bond_length * Sum_radii.\n"
      "By default the upper limit is set to 1.20, as a rule of thumb.\n"
      "The default should work for most crystalline materials, and covalent\n"
      "non-crystalline materials.\n"
      "Amorphous and liquids should increase this parameter to match the\n"
      "desired distance to cut_off the bonds\n\n"
      "-o, --out_file\n"
      "The output file name, by default the input seed name will be used.\n\n"
      "CREATOR:\n"
      "This program was created by PhD. Isaias Rodriguez Aguirre in May 2020.\n"
      "e-mail: isurwars@gmail.com\n"
      "This software was created during a Posdoctoral fellowship in IIM-UNAM.\n"
      "Thanks to DGAPA-UNAM for my financial support during my posdoctoral.\n"
      "And their continious support as part of the projects: IN104617 and \n"
      "IN116520.\n";
    std::cout << help_txt;

    exit(1);
} // PrintHelp

void ArgParser(int argc, char ** argv)
{
    const char * const short_opts = "r:b:o:w:h";
    const option long_opts[]      = {
        { "r_cut",       required_argument, nullptr, 'r' },
        { "bond_length", required_argument, nullptr, 'b' },
        { "out_file",    required_argument, nullptr, 'o' },
        { "bin_width",   required_argument, nullptr, 'w' },
        { "help",        no_argument,       nullptr, 'h' },
        { nullptr,       no_argument,       nullptr, 0   }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (opt == -1)
            break;

        switch (opt) {
            case 'r': // -r ot --r_cut
                _r_cut_ = std::stof(optarg);
                break;

            case 'b': // -b or --bond_length
                _bond_len_ = std::stof(optarg);
                break;

            case 'o': // -o or --out_file
                _out_file_name_ = optarg;
                break;

            case 'w': // -w or -- bin_width
                _bin_w_ = std::stof(optarg);
                break;

            case 'h': // -h or --help
                PrintHelp();
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
        std::cout << "Only one structure file most be provided." << '\n';
        exit(1);
    } else {
        if (argc > 1) {
            _in_file_name_ = argv[optind];
        } else {
            PrintHelp();
        }
    }
} // ProcessArgs

int main(int argc, char ** argv)
{
    Cell MyCell; // Struture to be analized
    std::list<Atom>::iterator MyAtom;
    std::pair<std::string, std::string> file_ext;

    /*
     * Parse the argumets provided as input
     */
    ArgParser(argc, argv);

    /*
     * Read the structure file
     */
    file_ext = GetExtension(_in_file_name_);
    if (_out_file_name_ == "") _out_file_name_ = file_ext.first;
    std::string MyExt = file_ext.second;
    // convert extension back to lower case
    std::for_each(MyExt.begin(), MyExt.end(), [](char & c){
        c = ::tolower(c);
    });
    if (MyExt == ".car") {
        MyCell = read_CAR(_in_file_name_);
    } else if (MyExt == ".cell") {
        MyCell = read_CELL(_in_file_name_);
    } else {
        std::cout << "File: " << MyExt << " currently not supported." << '\n';
        PrintHelp();
    }
    std::cout << "File " << _in_file_name_ << " opened successfully." << '\n';

    /*
     * This function calculate the distances from every atom to every
     * in the structure. A supercell method is used to create the images
     * in all directions up to _r_cut_ distance.
     */
    MyCell.RDF(_r_cut_, _bond_len_);

    /*
     * This function calculate the partial coordination number for pair of
     * elements. Bonded Atoms use the same parameters for BAD.
     */
    MyCell.CN();

    /*
     * This functions calculate the angle between every atom, and all pairs
     * of bonded atoms, the bonded atoms are calculated in Cell::RDF and it
     * must be called first.
     */
    MyCell.BAD();

    /*
     * This function use the distances to calculate J(r) and g(r).
     * The _r_cut_ parameters is the maximum distance,
     * the _bin_w_ parameter is the bin width to be used.
     */
    std::cout << "Writing output files: " << _out_file_name_ << '\n';
    MyCell.RDF_Histogram(_out_file_name_, _r_cut_, _bin_w_);
    MyCell.CN_Histogram(_out_file_name_);

    /*
     * This function use the angles to calculate the BAD.
     * The theta_max parameters is the maximum angle to compute,
     * the _bin_w_ parameter is the bin width to be used.
     */
    MyCell.BAD_Histogram(_out_file_name_);

    return 0;
} // main
