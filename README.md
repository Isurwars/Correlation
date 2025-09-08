![](Images/Banner.png)

# `Correlation`: An Analysis Tool for Liquids and for Amorphous Solids

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5514113.svg)](https://doi.org/10.5281/zenodo.5514113) [![Version](https://img.shields.io/badge/version-1.0.4-green)](https://img.shields.io/badge/version-1.0.4-green) [![License](https://img.shields.io/badge/license-MIT-brightgreen)](https://img.shields.io/badge/license-MIT-brightgreen) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02976/status.svg)](https://doi.org/10.21105/joss.02976)

`Correlation` is a high-performance, user-friendly tool for calculating and analyzing the structural properties of materials. It is designed for researchers working with atomistic simulations of liquids, amorphous solids, and crystalline structures.

The software computes key correlation functions from atomic coordinate files and exports the results in clean, ready-to-plot CSV files, making it easy to integrate into scientific workflows.

## Table of Contents
- [Key Features](#key-features)
- [Quick Start: Installation](#quick-start-installation)
  - [Prerequisites](#prerequisites)
  - [Windows](#windows)
  - [Linux (Debian/Ubuntu)](#linux-debian-ubuntu)
  - [Linux (Arch/Manjaro)](#linux-arch-manjaro)
  - [MacOS](#macos)
  - [Build Instructions](#build-instructions)
- [Usage](#usage)
- [Command-Line Options](#command-line-options)
- [License](#license)
- [Authors](#authors)
- [Acknowledgments](#acknowledgments)

## Key Features

Comprehensive Analysis: Calculates a full suite of standard correlation functions:

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

High Performance: The core calculation loops are parallelized using modern C++ techniques, enabling the analysis of systems with hundreds of thousands of atoms.

Data Smoothing: Includes built-in kernel smoothing (Gaussian, Triweight) to clean up noisy data for better presentation and analysis.

## Quick Start: Installation

### Prerequisites

- A modern C++ compiler (c++20 support required)
- CMake (version 3.20+)
- git
- Intel TBB (for parallelization)

### Windows

Installing MSYS2:
```
https://www.msys2.org/
```

Install prerequisites:

```bash
pacman -Syu
pacman -S --needed base-devel mingw-w64-x86_64-toolchain
pacman -S cmake tbb git
```

### Linux (Debian/Ubuntu):

```bash
sudo apt update
sudo apt install build-essential cmake tbb git
```

### Linux (Arch/Manjaro):

```bash
sudo pacman -Syu
sudo pacman base-devel cmake tbb git
```

### MacOS:

```bash
xcode-select --install
brew install cmake tbb git
```

### Build Instructions

#### Clone the repository:

```bash
git clone https://github.com/yourusername/correlation.git
cd correlation
```

#### Build the project:

```bash
rm -rf build && mkdir build && cd build
cmake ..
cmake --build .
```

#### Run tests:

```bash
ctest -V
```


#### (OPTIONAL) Install system-wide:

```bash
sudo cmake --install .
```



## Usage

The program is run from the command line, with the only required argument being the input structure file.

### Basic Command

```bash
./build/correlation <path_to_your_file>
```
### For example, to analyze a silicon structure:

```bash
./build/correlation silicon.car
```
This will run the analysis with default parameters and create output files (e.g., silicon_g.csv, silicon_PAD.csv) in the same directory as the input file.

### Example with Options

This will run the analysis with default parameters and create output files (e.g., silicon_g.csv, silicon_PAD.csv) in the same directory as the input file.
```bash
./build/correlation -R 10.0 -r 0.02 -S -K 0.05 -o si_run_1 si_crystal.car
```
This command will:

- Set the radial cutoff (-R) to 10.0 Å.
- Set the RDF bin width (-r) to 0.02 Å.
- Enable smoothing (-S).
- Set the smoothing kernel width (-K) to 0.05.
- Set the output file base name (-o) to si_run_1.

## Command-Line Options

| Option | Long Option | Argument | Description | Default |
|--------|-------------|----------|-------------|---------|
| -h | --help | - | Display the help text and exit. | - |
| -o | --out_file | <path> | The base name for output files.| Input filename |
| -R | --r_cut | <float> | Cutoff radius for RDF calculations (in Angstroms). | 20.0 |
| -r | --r_bin_width | <float> | Width of the histogram bins for RDFs (in Angstroms). | 0.05 |
| -a | --angle_bin_width | <float> | Width of the histogram bins for PAD (in degrees). | 1.0 |
| -S | --smoothing | - | Enable kernel smoothing on all calculated distributions. | Disabled |
| -k | --kernel | <int> | Smoothing kernel type: 1=Gaussian, 2=Bump, 3=Triweight. | 1 (Gaussian) |
| -K | --kernel_sigma | <float> | Width (sigma) of the smoothing kernel. | 0.081 |


## Built with

- [emacs](https://www.gnu.org/software/emacs/) - An extensible, customizable, free/libre text editor — and more.
- [MSYS2](https://www.msys2.org/) - Software Distribution and Building Platform for Windows

## Authors

- **Isaías Rodríguez** - _Corresponding Author_ - [Isurwars](https://github.com/Isurwars) <isurwars@gmail.com>
- **Salvador Villareal Lopez Rebuelta** <salvadorvillarreallr@gmail.com>
- **Renela M. Valladares** <renelavalladares@gmail.com>
- **Alexander Valladares** <valladar@ciencias.unam.mx>
- **David Hinojosa-Romero** <david18_hr@ciencias.unam.mx>
- **Ulises Santiago**
- **Ariel A. Valladares** <valladar@unam.mx>

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

I.R. acknowledge PAPIIT, DGAPA-UNAM for his postdoctoral fellowship.
D.H.R. acknowledge Consejo Nacional de Ciencia y Tecnología (CONACyT) for supporting his graduate studies.
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under Grant No. IN104617 and IN116520.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out in the Supercomputing Center of DGTIC-UNAM.
I.R. would like to express his gratitude to F. B. Quiroga, M. A. Carrillo, R. S. Vilchis, S. Villareal and A. de Leon, for their time invested in testing the code, as well as the structures provided for benchmarks and tests.
