![](Images/Banner.png)

# `Correlation`: An Analysis Tool for Liquids and for Amorphous Solids

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5514113.svg)](https://doi.org/10.5281/zenodo.5514113) [![Version](https://img.shields.io/badge/version-1.0.4-green)](https://img.shields.io/badge/version-1.0.4-green) [![License](https://img.shields.io/badge/license-MIT-brightgreen)](https://img.shields.io/badge/license-MIT-brightgreen) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02976/status.svg)](https://doi.org/10.21105/joss.02976)

`Correlation` is an analysis tool for correlation functions and correlation related properties of materials. In particular, for atomistic structure files of heavily used material simulation software like: DMoL3 (_.CAR), CASTEP(_.CELL), ONETEP(_.DAT), LAMMPS(_.XYZ),etc...

## Table of Contents
- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Build Instructions](#build-instructions)
- [Usage](#usage)
- [License](#license)
- [Authors](#authors)
- [Acknowledgments](#acknowledgments)

## Features

This program calculates the main correlation functions of a material:

- Radial Distribution Function (J(r))
- Pair Distribution Function (g(r))
- Reduced Pair Distribution Function (G(r))
- Coordination Number (CN)
- Plane-Angle Distribution (PAD)
- Structure Factor (S(Q))

Supports structure files from:
- DMoL3 (`.CAR`)
- CASTEP (`.CELL`)
- ONETEP (`.DAT`)
- LAMMPS (`.XYZ`)


## Installation

### Prerequisites

#### Windows

There are several (GCC) implementations for Windows.
WE recommend the MSYS2 implementation:

Installing MSYS2:
```
MSYS2: https://www.msys2.org/
```

Using MSYS2 (recommended):
```
pacman -Syu
pacman -S --needed base-devel mingw-w64-x86_64-toolchain
pacman -S cmake tbb git
```

#### Linux
##### Debian/Ubuntu:

```
sudo apt update
sudo apt install build-essential cmake tbb git
```

##### Arch/Manjaro:
```
sudo pacman -Syu
sudo pacman base-devel cmake tbb git
```

#### MacOS
```
xcode-select --install
brew install cmake tbb git
```

### Building from Source
```
# Clean build and compile
rm -rf build && mkdir build && cd build
cmake ..
cmake --build .

# Install system-wide (optional)
sudo cmake --install .
```



## Usage

---

```
    USAGE: correlation [OPTIONS] [input_file]

      The minimal argument is a structure file, this program requires a file
      that contains atom positions, crystal structure and composition.
      Supported structure files are:
        -*.CAR   Materials Studio structure file.
        -*.CELL  CASTEP structure file.
        -*.DAT   ONETEP structure file.
        -*.XYZ   LAMMPS structure file.

      OPTIONS:
        HELP OPTIONS
          -h, --help
            Display this help text.

        RADIAL OPTIONS:
          -n, --normalize
            Used to switch between weighted partials (default), or normalize all the partials to 1 when r tends to infinity.
          -r, --r_cut
            Cutoff radius in the calculation of g(r), G(r) and J(r). The default
            radius it's set to 2 nm. The maximum recommended radius is the same as
            the shortest length of the lattice parameters of the cell, anything
            above this PBC value can be affected by periodic interactions.
          -s, --self_interaction
            Include self-interactions, by default false.
          -w, --bin_width
            Width of the histograms for g(r) and J(r), the default is 0.05 nm.

        STRUCTURE FACTOR OPTIONS:
          -q, --q_bin_width
            Width of the histograms for S(Q), the default is 0.157079 nm^()-1.)

        BOND-ANGLE OPTIONS:
          -a, --angle_bin_width
            Width of the histograms for the PAD, default set to 1.0°.
          -b, --bond_parameter
            The ideal covalent bond length is the sum of covalent radii
            of the two atoms. The criterion used to consider atoms as bonded
            is the following:
              0.6 * Sum_radii < distance < bond_parameter * Sum_radii.
            By default the bond_parameter is set to 1.30, as a rule of thumb.
            The default should work for most crystalline materials,
            as well as most covalent non-crystalline materials.
            For amorphous and liquid materials the bond_parameter should be
            increased to match the desired distance to cut_off the bonds.
            A bond_parameter of 1.42 is recomended for amorphous materials.
          -i, --in_bond_file
            The input file with the bond distances for every pair of elements
            in the corresponding input structure. The file should have the
            following format:
              Si Si 2.29
              Mg Mg 2.85
              C  C  1.55
              C  Si 1.86
              Si Mg 2.57
              C  Mg 2.07
            If any of the pairs is missing in the input file, the corresponding
            bond distance will be set using the bond_parameter(1.30 by default).

        OUTPUT OPTIONS:
          -o, --out_file
            The output file name, by default the input seed name will be used.

        SMOOTHING OPTIONS:
          -k, --kernel
            Smoothing kernel selector (Default: 1):
              1: Gaussian kernel.
              2: Bump Function kernel.
              3: Triweight kernel.
          -K, --kernel_sigma
            Width of the smoothing kernel, by default 0.081.
          -S, --smoothing
            Smoothing is disabled by default, this option enable smoothing.



## Built With

- [emacs](https://www.gnu.org/software/emacs/) - An extensible, customizable, free/libre text editor — and more.
- [MSYS2](https://www.msys2.org/) - Software Distribution and Building Platform for Windows

## Authors

- **Isaías Rodríguez** - _Corresponding Author_ - [Isurwars](https://github.com/Isurwars) <isurwars@gmail.com>
- **Renela M. Valladares** <renelavalladares@gmail.com>
- **Alexander Valladares** <valladar@ciencias.unam.mx>
- **David Hinojosa-Romero** <david18_hr@ciencias.unam.mx>
- **Ulises Santiago**
- **Ariel A. Valladares** <valladar@unam.mx>

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details


- Support for additional output files, like hdf5 standard.

- Inclusion of other correlation functions like Velocity Correlation Functions, to further improve the analysis of liquids and phase transitions.


# Acknowledgments

I.R. acknowledge PAPIIT, DGAPA-UNAM for his postdoctoral fellowship.
D.H.R. acknowledge Consejo Nacional de Ciencia y Tecnología (CONACyT) for supporting his graduate studies.
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under Grant No. IN104617 and IN116520.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out in the Supercomputing Center of DGTIC-UNAM.
I.R. would like to express his gratitude to F. B. Quiroga, M. A. Carrillo, R. S. Vilchis, S. Villareal and A. de Leon, for their time invested in testing the code, as well as the structures provided for benchmarks and tests.
