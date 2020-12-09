# **Correlation**: An Analyzing Tool for Liquids and for Amorphous Solids

**Correlation** is an analyzing tool of correlation functions and correlation related properties of materials. In particular, for atomistic structure files of heavily used material simulation software like: DMoL3 (.CAR), CASTEP(.CELL), ONETEP(*.DAT), LAAMPS(.XYZ),etc...  

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This projects is designed to be compatible with GNU Compiler Collection (GCC).
```
Windows:

There are several (GCC) implementations for Windows, developers recommend the following implementations:

MinGW: http://www.mingw.org/

MSYS2: https://www.msys2.org/

Pre-Compiled binary for Windows 10 is present in the /bin/ folder.
```

```
Unix:

The code was tested with Ubuntu and Debian, but it should work on other distros.
```


```
MacOS:

MacOS is currently not supported. (Major priority for 1.0 release)
Clang , the default c/c++ compiler,  doesn't compile even with some c++17 flags activated.
```



### Installing

For compiling you should have GCC environment correctly installed on the system.

To compile the program:

```
make all
```

And then to clean the precompiled objects files:

```
make clean
```

The newly compiled program is located in the bin folder.

## Running the tests

Three Atomistic Structure Files are included for testing.

You can run all three tests with the command:

```
make tests
```

### Break down into end to end tests

The first test include one of the structure file for amorphous palladium as published by Rodríguez et al. https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.024422

To run this test you could run:
```
make test1
```
or

```
./bin/correlarion.exe ./test/test_1/aPd.cell
```

Three result(s) file(s) should be generated in the same directory as the input fil(aPd.cell).
These Coma Separeted Values files can then be analyzed with your favorite tool like: LibreOffice Math, Office Excel, OriginPro, etc...

## Built With

* [Atom](https://atom.io/) - A hackable text editor for the 21st Century
* [MSYS2](https://www.msys2.org/) - Software Distribution and Building Platform for Windows


## Authors

* **Isaías Rodríguez** - *Initial work* - [Isurwars](https://github.com/Isurwars)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

# Acknowledgments

I.R. acknowledge PAPIIT, DGAPA-UNAM for his posdoctoral fellowship.
D.H.R. acknowledge Consejo Nacional de Ciencia y Tecnología (CONACyT) for supporting his graduate studies.
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under Grant No. IN104617 and IN116520.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out in the Supercomputing Center of DGTIC-UNAM.
I.R. would like to express his gratitude to F. B. Quiroga, M. A. Carrillo, R. S. Vilchis, S. Villareal and A. de Leon, for their time invested in testing the code, as well as the structures provided for benchmarks and tests.
