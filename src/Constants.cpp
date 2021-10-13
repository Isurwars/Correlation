/*
 * Copyright [2021] <@isurwars>
 */
#include "Constants.h"

/*
 * We define a simple function that returns the covalent radii.
 * Covalent Radii definitions in Angstroms.
 * Based on the following sources:
 * Molecular single-bond Covalent Radii for elements 1-118.
 * Pyykko, P. and Atsumi, M., DOI: 10.1002/chem.200800987
 * Colvant Radii revisited, Cordero B, et al.
 * DOI: 10.1039/b801115j
 */

double Covalent_Radii(std::string element) {
  /*1s*/
  if (element == "H") return 0.32;

  if (element == "He") return 0.46;

  /*2s*/
  if (element == "Li") return 1.33;

  if (element == "Be") return 1.02;

  /*2p*/
  if (element == "B") return 0.85;

  if (element == "C") return 0.75;

  if (element == "N") return 0.71;

  if (element == "O") return 0.63;

  if (element == "F") return 0.64;

  if (element == "Ne") return 0.67;

  /*3s*/
  if (element == "Na") return 1.55;

  if (element == "Mg") return 1.39;

  /*3p*/
  if (element == "Al") return 1.26;

  if (element == "Si") return 1.16;

  if (element == "P") return 1.11;

  if (element == "S") return 1.03;

  if (element == "Cl") return 0.99;

  if (element == "Ar") return 0.96;

  /*4s*/
  if (element == "K") return 1.96;

  if (element == "Ca") return 1.71;

  /*3d*/
  if (element == "Sc") return 1.48;

  if (element == "Ti") return 1.36;

  if (element == "V") return 1.34;

  if (element == "Cr") return 1.22;

  if (element == "Mn") return 1.19;

  if (element == "Fe") return 1.16;

  if (element == "Co") return 1.11;

  if (element == "Ni") return 1.10;

  if (element == "Cu") return 1.12;

  if (element == "Zn") return 1.18;

  /*4p*/
  if (element == "Ga") return 1.24;

  if (element == "Ge") return 1.21;

  if (element == "As") return 1.21;

  if (element == "Se") return 1.16;

  if (element == "Br") return 1.14;

  if (element == "Kr") return 1.17;

  /*5s*/
  if (element == "Rb") return 2.10;

  if (element == "Sr") return 1.85;

  /*4d*/
  if (element == "Y") return 1.63;

  if (element == "Zr") return 1.54;

  if (element == "Nb") return 1.47;

  if (element == "Mo") return 1.38;

  if (element == "Tc") return 1.28;

  if (element == "Ru") return 1.25;

  if (element == "Rh") return 1.25;

  if (element == "Pd") return 1.20;

  if (element == "Ag") return 1.28;

  if (element == "Cd") return 1.36;

  /*5p*/
  if (element == "In") return 1.42;

  if (element == "Sn") return 1.40;

  if (element == "Sb") return 1.40;

  if (element == "Te") return 1.36;

  if (element == "I") return 1.33;

  if (element == "Xe") return 1.31;

  /*6s*/
  if (element == "Cs") return 2.32;

  if (element == "Ba") return 1.96;

  /*4f*/
  if (element == "Ce") return 1.63;

  if (element == "Pr") return 1.76;

  if (element == "Nd") return 1.74;

  if (element == "Pm") return 1.73;

  if (element == "Sm") return 1.72;

  if (element == "Eu") return 1.68;

  if (element == "Gd") return 1.69;

  if (element == "Tb") return 1.68;

  if (element == "Dy") return 1.67;

  if (element == "Ho") return 1.66;

  if (element == "Er") return 1.65;

  if (element == "Tm") return 1.64;

  if (element == "Yb") return 1.70;

  if (element == "Lu") return 1.62;

  /*5d*/
  if (element == "La") return 1.80;

  if (element == "Hf") return 1.52;

  if (element == "Ta") return 1.46;

  if (element == "W") return 1.37;

  if (element == "Re") return 1.31;

  if (element == "Os") return 1.29;

  if (element == "Ir") return 1.22;

  if (element == "Pt") return 1.23;

  if (element == "Au") return 1.24;

  if (element == "Hg") return 1.33;

  /*6p*/
  if (element == "Tl") return 1.44;

  if (element == "Pb") return 1.44;

  if (element == "Bi") return 1.51;

  if (element == "Po") return 1.45;

  if (element == "At") return 1.47;

  if (element == "Rn") return 1.42;

  /*7s*/
  if (element == "Fr") return 2.23;

  if (element == "Ra") return 2.01;

  /*5f*/
  if (element == "Th") return 1.75;

  if (element == "Pa") return 1.69;

  if (element == "U") return 1.70;

  if (element == "Np") return 1.71;

  if (element == "Pu") return 1.72;

  if (element == "Am") return 1.66;

  if (element == "Cm") return 1.66;

  if (element == "Bk") return 1.68;

  if (element == "Cf") return 1.68;

  if (element == "Es") return 1.65;

  if (element == "Fm") return 1.67;

  if (element == "Md") return 1.73;

  if (element == "No") return 1.76;

  if (element == "Lr") return 1.61;

  /*6d*/
  if (element == "Ac") return 1.86;

  if (element == "Rf") return 1.57;

  if (element == "Db") return 1.49;

  if (element == "Sg") return 1.43;

  if (element == "Bh") return 1.41;

  if (element == "Hs") return 1.34;

  if (element == "Mt") return 1.29;

  if (element == "Ds") return 1.28;

  if (element == "Rg") return 1.21;

  if (element == "Cn") return 1.22;

  /*7p*/
  if (element == "Nh") return 1.36;

  if (element == "Fl") return 1.43;

  if (element == "Mc") return 1.62;

  if (element == "Lv") return 1.75;

  if (element == "Ts") return 1.65;

  if (element == "Og") return 1.57;

  return 0;
}  // Covalent_Radii
