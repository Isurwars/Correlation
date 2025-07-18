![](Images/Banner.png)

# `Correlation`: An Analysis Tool for Liquids and for Amorphous Solids

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5514113.svg)](https://doi.org/10.5281/zenodo.5514113) [![Version](https://img.shields.io/badge/version-1.0.4-green)](https://img.shields.io/badge/version-1.0.4-green) [![License](https://img.shields.io/badge/license-MIT-brightgreen)](https://img.shields.io/badge/license-MIT-brightgreen) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02976/status.svg)](https://doi.org/10.21105/joss.02976)

`Correlation` is an analysis tool for correlation functions and correlation related properties of materials. In particular, for atomistic structure files of heavily used material simulation software like: DMoL3 (_.CAR), CASTEP(_.CELL), ONETEP(_.DAT), LAMMPS(_.XYZ),etc...

## Table of Contents
- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Build Instructions](#build)
- [Usage](#usage)
- [License](#license)
- [Authors](#authors)
- [Acknowledgments](#acknowledgments)

## Features

Correlation calculates essential correlation functions for material analysis:

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

#### Windows (MSYS2 recommended)

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

#### Linux (Debian/Ubuntu):

```bash
sudo apt update
sudo apt install build-essential cmake tbb git
```

#### Linux (Arch/Manjaro):

```bash
sudo pacman -Syu
sudo pacman base-devel cmake tbb git
```

#### MacOS:

```bash
xcode-select --install
brew install cmake tbb git
```

### Build

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

---
---


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
