![](Images/Banner.png)

# `Correlation`: An Analysis Tool for Liquids and for Amorphous Solids

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5514113.svg)](https://doi.org/10.5281/zenodo.5514113) [![Version](https://img.shields.io/badge/version-1.8.2-green)](https://img.shields.io/badge/version-1.8.2-green) [![License](https://img.shields.io/badge/license-MIT-brightgreen)](https://img.shields.io/badge/license-MIT-brightgreen) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02976/status.svg)](https://doi.org/10.21105/joss.02976)

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
- [License](#license)
- [Authors](#authors)
- [Acknowledgments](#acknowledgments)

## Key Features

Comprehensive Analysis: Calculates a full suite of standard correlation
functions:

- Radial Distribution Function (J(r))
- Pair Distribution Function (g(r))
- Reduced Pair Distribution Function (G(r))
- Coordination Number (CN)
- Plane-Angle Distribution (PAD)
- Structure Factor (S(Q))
- Velocity Autocorrelation Function (VACF)
- Vibrational Density of States (VDOS)


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
- Slint (for GUI)

### Windows

1.  **Install Visual Studio**: Download and install [Visual Studio](https://visualstudio.microsoft.com/downloads/) with the "Desktop development with C++" workload. This includes the MSVC compiler and CMake.
2.  **Install Rust**: Download and install Rust from [rust-lang.org](https://www.rust-lang.org/tools/install).
3.  **Install Git**: Download and install Git from [git-scm.com](https://git-scm.com/download/win).

### Linux (Debian/Ubuntu)

```bash
sudo apt update
sudo apt install build-essential cmake git rustc cargo
# Optional: Install TBB and HDF5 if you prefer system libs over fetching
sudo apt install libtbb-dev libhdf5-dev
```

### Linux (Arch/Manjaro)

```bash
sudo pacman -Syu
sudo pacman -S base-devel cmake git rust
# Optional: Install TBB and HDF5 if you prefer system libs over fetching
sudo pacman -S intel-tbb hdf5
```

### MacOS

```bash
xcode-select --install
brew install cmake git rustup
# Optional: Install TBB and HDF5 if you prefer system libs over fetching
brew install tbb hdf5
```

### Build Instructions

#### Clone the repository:

```bash
git clone https://github.com/isurwars/correlation.git
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

`Correlation` features an intuitive graphical user interface (GUI) to guide you through the analysis process. To start the application, run the executable:

```bash
./build/correlation
```

### 1. Load a Structure File
Click the **"Load a structure file"** button in the **Input File** card. This opens a file dialog to select your material structure. Supported formats include `.car`, `.cell`, `.dat`, and `.xyz`.

### 2. File Information
Once a file is loaded, the **File Info** card displays the total atom count and the breakdown by element type.

### 3. Configure Analysis Options
Adjust the calculation parameters in the **Options** card:
- **RDF Max r & Bin Width:** Set the maximum radius ($r_{max}$) and the bin width for Radial Distribution Functions.
- **Max Q & S(Q) Bin Width:** Set the maximum scattering vector ($Q_{max}$) and the bin width for the Structure Factor.
- **FFT Max r:** Set the maximum integration radius for the Fourier Transform.
- **PAD Bin Width:** Set the bin width for the Plane-Angle Distribution.
- **Smoothing:** Enable kernel smoothing and select the kernel type (Gaussian, Bump, or Triweight) and sigma value.

### 4. Bond Cutoffs
The **Bond Cutoffs** card allows you to review and manually adjust the distances used to define atomic bonds, which are used for coordination number and angle calculations.

### 5. Run Analysis
Configure the final output and execution settings in the **Run Analysis** card:
-   **Export Format**: Choose to export results as **CSV** (for easy plotting) or **HDF5** (for large datasets/archival).
-   **Frame Selection**: Specify the **Start Frame** and **End Frame** to analyze a specific subset of your trajectory.
-   **Run Analysis**: Click the **"Run Analysis"** button to start the computations. Progress will be displayed in the bar below.
-   **Write Files**: Once the analysis is complete, click **"Write Files"** to save the results to disk.
## Built with

- [emacs](https://www.gnu.org/software/emacs/) - An extensible, customizable, free/libre text editor — and more.

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
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under Grant No. IN104617, IN116520 and IIN118223.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out in the Supercomputing Center of DGTIC-UNAM.
I.R. would like to express his gratitude to F. B. Quiroga, M. A. Carrillo, R. S. Vilchis, S. Villareal and A. de Leon, for their time invested in testing the code, as well as the structures provided for benchmarks and tests.
