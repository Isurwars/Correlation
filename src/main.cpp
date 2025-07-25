// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <getopt.h> // For the argument parsing

#include <algorithm>  // for_each
#include <filesystem> // Getting extension, making subdirectories
#include <iostream>   // Standar IO library
#include <string>     // String manipulation and algorithms

#include "../include/Cell.hpp"       // Cell Class
#include "../include/ReadFiles.hpp"  // File reading and parsing
#include "../include/WriteFiles.hpp" // Write Output files

inline std::string _in_file_name_ = "";
inline std::string _out_file_name_ = "";
inline std::string _bond_file_name_ = "";
bool _bond_in_file_ = false;
bool _normalize_ = false;
bool _self_interaction_ = false;
bool _smoothing_ = false;
double _r_cut_ = 20.0;
double _bin_w_ = 0.050000;
double _q_bin_w_ = 0.15707963267949;
double _q_cut_ = 180.0;
double _bond_par_ = 1.3;
double _angle_bin_w_ = 1.00000;
double _smooth_sigma_ = 0.081;
int _kernel_ = 1;

void PrintHelp() {
  const std::string help_txt =
      "CORRELATION\n"
      "DESCRIPTION:\n"
      "  This program calculates the main correlation functions of a "
      "material:\n"
      "    - Radial Distribution Function (J(r)).\n"
      "    - Pair Distribution Function (g(r)).\n"
      "    - Reduced Pair Distribution Function (G(r)).\n"
      "    - Coordination Number (CN).\n"
      "    - Plane-Angle Distribution (PAD).\n\n"
      "    - Structure Factor (S(Q)).\n"
      "USAGE: correlation [OPTIONS] [input_file]\n"
      "  The minimal argument is a structure file, this program requires a "
      "file\n"
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
      "      -n, --normalize\n"
      "        Used to switch between weighted partials (default), or "
      "normalize\n"
      "        all the partials to 1 when r tends to infinity\n"
      "      -R, --r_cut\n"
      "        Cutoff radius in the calculation of g(r), G(r) and J(r). The "
      "default\n"
      "        radius it's set to 2 nm. The maximum recommended radius is the "
      "same as\n"
      "        shortest length of the periodic boundary conditions (PBC), "
      "anything\n"
      "        above this PBC value can be affected by periodic "
      "interactions.\n\n"
      "      -s, --self_interaction\n"
      "        Include self-interactions, by default false.\n"
      "      -r, --r_bin_width\n"
      "        Width of the histograms for g(r) and J(r), the default is 0.05 "
      "nm.\n\n"
      "    STRUCTURE FACTOR OPTIONS:\n"
      "      -q, --q_bin_width\n"
      "        Width of the histograms for S(Q), the default is 0.157079 "
      "nm^(-1).\n\n"
      "      -Q, --q_cut\n "
      "        Max vector Q for S(Q), the default is 180.0.\n\n"
      "    BOND-ANGLE OPTIONS:\n"
      "      -a, --angle_bin_width\n"
      "        Width of the histograms for the PAD, default set to 1.0°.\n\n"
      "      -b, --bond_parameter\n"
      "        The ideal covalent bond length is the sum of covalent radii\n"
      "        of the two atoms. The criterion used to consider atoms as "
      "bonded\n"
      "        is the following:\n"
      "            0.6 * Sum_radii < distance < bond_parameter * Sum_radii.\n"
      "        By default the bond_parameter is set to 1.30, as a rule of "
      "thumb.\n"
      "        The default should work for most crystalline materials,\n"
      "        as well as most covalent non-crystalline materials.\n"
      "        For amorphous and liquid materials the bond_parameter should "
      "be\n"
      "        increased to match the desired distance to cut_off the bonds. \n"
      "        A bond_parameter of 1.42 is recomended for amorphous "
      "materials.\n\n"
      "      -i, --in_bond_file"
      "        The input file with the bond distances for every pair of "
      "elements\n"
      "        in the corresponding input structure. The file should have the\n"
      "        following format:\n"
      "             Si Si 2.29\n"
      "             Mg Mg 2.85\n"
      "             C  C  1.55\n"
      "             C  Si 1.86\n"
      "             Si Mg 2.57\n"
      "             C  Mg 2.07\n"
      "        If any of the pairs is missing in the input file, the "
      "corresponding\n"
      "        bond distance will be set using the bond_parameter(1.30 by "
      "default).\n\n"
      "    OUTPUT OPTIONS:\n"
      "      -o, --out_file\n"
      "        The output file name, by default the input seed name will be "
      "used.\n\n"
      "    SMOOTHING OPTIONS:\n"
      "      -k, --kernel\n"
      "        Smoothing kernel selector (Default: 1):\n"
      "          1: Gaussian kernel.\n"
      "          2: Bump Function kernel.\n"
      "          3: Triweight kernel.\n"
      "      -K, --kernel_sigma\n"
      "        Width of the smoothing kernel, by default 0.081.\n"
      "        The recommended values for sigma are:\n"
      "        Gaussian Kernel: 0.081\n"
      "        Bump Kernel: 0.19\n"
      "      -S, --smoothing\n"
      "        Smoothing is disabled by default, this option enable "
      "smoothing\n\n"
      "CREATOR:\n"
      "  This program was created by PhD. Isaias Rodriguez Aguirre, November "
      "2020.\n"
      "  e-mail: isurwars@gmail.com\n"
      "ACKNOWLEDGMENTS:\n"
      "  This software was created during a Posdoctoral fellowship in "
      "IIM-UNAM.\n"
      "  Thanks to DGAPA-UNAM for the financial support during my fellowship.\n"
      "  And their continious support as part of the Workgroup projects with "
      "IDs:\n"
      "  IN104617 and IN116520.\n";
  std::cout << help_txt;

  exit(1);
} // PrintHelp

void ArgParser(int argc, char **argv) {
  const char *const short_opts = "a:b:hi:k:K:no:q:Q:r:R:sS";
  const option long_opts[] = {
      {"angle_bin_width", required_argument, nullptr, 'a'},
      {"bond_parameter", required_argument, nullptr, 'b'},
      {"help", no_argument, nullptr, 'h'},
      {"in_bond_file", required_argument, nullptr, 'i'},
      {"kernel", required_argument, nullptr, 'k'},
      {"kernel_sigma", required_argument, nullptr, 'K'},
      {"normalize", no_argument, nullptr, 'n'},
      {"out_file", required_argument, nullptr, 'o'},
      {"q_bin_width", required_argument, nullptr, 'q'},
      {"q_cut", required_argument, nullptr, 'Q'},
      {"r_bin_width", required_argument, nullptr, 'r'},
      {"r_cut", required_argument, nullptr, 'R'},
      {"self_interaction", no_argument, nullptr, 's'},
      {"smoothing", no_argument, nullptr, 'S'},
      {nullptr, no_argument, nullptr, 0}};

  while (true) {
    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (opt == -1)
      break;

    switch (opt) {
    case 'a': // -a or --angle_bin_width
      try {
        _angle_bin_w_ = std::stof(optarg);
      } catch (const std::exception &e) {
        std::cout << "Invalid input argument: '" << optarg
                  << "' in angle_bin_width parameter '-a' "
                  << "(Real number expected)." << std::endl;
        exit(1);
      }
      break;
    case 'b': // -b or --bond_parameter
      try {
        _bond_par_ = std::stof(optarg);
      } catch (const std::exception &e) {
        std::cout << "Invalid input argument: '" << optarg
                  << "' in bond_parameter '-b' "
                  << "(Real number expected)." << std::endl;
        exit(1);
      }
      break;
    case 'h': // -h or --help
      PrintHelp();
      break;
    case 'i': // -i or --in_bond_file
      _bond_in_file_ = true;
      _bond_file_name_ = optarg;
      break;
    case 'k': // -k or --kernel
      _kernel_ = std::stoi(optarg);
      break;
    case 'K': // -K or --kernel_sigma
      _smooth_sigma_ = std::stof(optarg);
      break;
    case 'n': // -n or --normalize
      _normalize_ = true;
      break;
    case 'o': // -o or --out_file
      _out_file_name_ = optarg;
      break;
    case 'q': // -q or --q_bin_width
      try {
        _q_bin_w_ = std::stof(optarg);
      } catch (const std::exception &e) {
        std::cout << "Invalid input argument: '" << optarg
                  << "' in q_bin_width parameter '-q' "
                  << "(Real number expected)." << std::endl;
        exit(1);
      }
      break;
    case 'Q': // -Q or --q_cut
      try {
        _q_cut_ = std::stof(optarg);
      } catch (const std::exception &e) {
        std::cout << "Invalid input argument: '" << optarg
                  << "' in q_cut parameter '-Q' "
                  << "(Real number expected)." << std::endl;
        exit(1);
      }
      break;
    case 'r': // -r or --r_bin_width
      try {
        _bin_w_ = std::stof(optarg);
      } catch (const std::exception &e) {
        std::cout << "Invalid input argument: '" << optarg
                  << "' in bin_width parameter '-w' "
                  << "(Real number expected)." << std::endl;
        exit(1);
      }
      break;
    case 'R': // -R or --r_cut
      try {
        _r_cut_ = std::stof(optarg);
      } catch (const std::exception &e) {
        std::cout << "Invalid input argument: '" << optarg
                  << "' in r_cut parameter '-r' "
                  << "(Real number expected)." << std::endl;
        exit(1);
      }
      break;
    case 's': // -s or --self_interaction
      _self_interaction_ = true;
      break;
    case 'S': // -S or --smoothing
      _smoothing_ = true;
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

int main(int argc, char **argv) {
  Cell MyCell; // Struture to be analized
  // std::list<Atom>::iterator MyAtom;
  // std::pair<std::string, std::string> file_ext;

  /*
   * Parse the input argumets provided
   */
  ArgParser(argc, argv);

  /*
   * Read the structure file
   */
  std::filesystem::path _in_file_{_in_file_name_};
  std::filesystem::path _out_file_;
  std::filesystem::path _in_path_ = _in_file_.parent_path();

  if (_out_file_name_.empty()) {
    _out_file_ = _in_path_ / _in_file_.stem();
    _out_file_name_ = _out_file_.generic_string();
  }

  std::string MyExt = _in_file_.extension().generic_string();
  // Convert extension back to lower case
  std::for_each(MyExt.begin(), MyExt.end(), [](char &c) { c = ::tolower(c); });
  if (MyExt == ".car") {
    std::cout << "Reading CAR file: " << _in_file_name_ << std::endl;
    MyCell = readCar(_in_file_name_);
  } else if (MyExt == ".cell") {
    std::cout << "Reading CELL file: " << _in_file_name_ << std::endl;
    MyCell = readCell(_in_file_name_);
  } else if (MyExt == ".dat") {
    std::cout << "Reading DAT file: " << _in_file_name_ << std::endl;
    MyCell = readOnetepDat(_in_file_name_);
  } else {
    std::cout << "File: " << MyExt << " currently not supported." << std::endl;
    PrintHelp();
  }
  std::cout << "File " << _in_file_name_ << " opened successfully."
            << std::endl;
  // Create Bond distance Matrix
  MyCell.populateBondLength(_bond_par_);
  if (_bond_in_file_) {
    /*Read the Bond Distances from the external file*/
    MyCell.setBondLength(readBond(_bond_file_name_, MyCell));
  }
  std::cout << "Bond lenght populated correctly" << std::endl;
  MyCell.populateElementID();
  std::cout << "Element IDs populated correctly" << std::endl;
  /*
   * This function calculates the distances between every pair of atoms
   * in the structure. A supercell method is used to create the images
   * in all directions up to _r_cut_ distance.
   */
  MyCell.distancePopulation(_r_cut_, _self_interaction_);
  std::cout << "Distance tensor populated correctly" << std::endl;
  /*
   * This function calculates the angle between every atom and all pairs
   * of bonded atoms. The bonded atoms are calculated in Cell::RDF and it
   * must be called first.
   */
  MyCell.planeAnglePopulation();
  std::cout << "Angle tensorpopulated correctly" << std::endl;

  DistributionFunctions MyDF(MyCell);
  /*
   * This function calculates the partial coordination number for pairs of
   * elements. Bonded Atoms use the same parameters for PAD.
   */
  MyDF.coordinationNumber();

  /*
   * This function uses the distances to calculate g(r), G(r) and J(r).
   * The _r_cut_ parameter is the cutoff distance,
   * the _bin_w_ parameter is the bin width to be used.
   */
  std::cout << "Calculating Coordination Functions: " << _out_file_name_
            << std::endl;
  MyDF.calculateRDF(_r_cut_, _bin_w_, _normalize_);

  /*
   * This function uses the angles to calculate the PAD.
   * The theta_max parameter is the maximum angle to compute,
   * the _bin_w_ parameter is the bin width to be used.
   */
  MyDF.calculatePAD(180.0, _angle_bin_w_);

  /*
   * This function calculate S(Q) as the Fourier Transform of G(r).
   * The _r_cut_ parameter is the cutoff distance in r-space,
   * the _q_bin_w_ parameter is the bin width to be used in q-space.
   */
  MyDF.calculateSQ(_q_cut_, _q_bin_w_, _normalize_);
  // MyCell.calculateSQExact(_q_cut_, _q_bin_w_, _normalize_);

  /*
   * This function calculate calculateXRD with bragg equation.
   * The lambda parameter is the wavelenght of the X-ray, default is 1.5406,
   * corresponding to Cu K_alpha.
   * The theta_min and theta_max parameters are the range to be calculated.
   * The bin_width paramater is the bin width to be used in theta.
   *
   * WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!
   *
   * This function is currently being tested and results should be taken
   * with a grain of salt!
   *
   * WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!
   *
   */
  MyDF.calculateXRD(1.5406, 10.0, 90.0, 0.2);

  if (_smoothing_) {
    std::cout << "Smoothing... " << std::endl;
    MyDF.Smoothing(_smooth_sigma_, _kernel_);
  }

  std::cout << "Writing output files: " << _out_file_name_ << std::endl;
  WriteCSV(MyDF, _out_file_name_, _smoothing_);

  std::cout << "Job in " << _in_file_name_ << " finished successfully." << '\n';
  return 0;

  std::vector<double> aux(MyDF.g().size(), 0);
} // main
