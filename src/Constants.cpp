/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2024 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include "Constants.hpp"

/*
 * We define a simple function that returns the covalent radii.
 * Covalent Radii definitions in Angstroms.
 * Based on the following sources:
 * Molecular single-bond Covalent Radii for elements 1-118.
 * Pyykko, P. and Atsumi, M., DOI: 10.1002/chem.200800987
 * Colvant Radii revisited, Cordero B, et al.
 * DOI: 10.1039/b801115j
 */

double covalentRadii(std::string_view element) {
  /*1s*/
  if (element == "H")
    return 0.32;

  if (element == "He")
    return 0.46;

  /*2s*/
  if (element == "Li")
    return 1.33;

  if (element == "Be")
    return 1.02;

  /*2p*/
  if (element == "B")
    return 0.85;

  if (element == "C")
    return 0.75;

  if (element == "N")
    return 0.71;

  if (element == "O")
    return 0.63;

  if (element == "F")
    return 0.64;

  if (element == "Ne")
    return 0.67;

  /*3s*/
  if (element == "Na")
    return 1.55;

  if (element == "Mg")
    return 1.39;

  /*3p*/
  if (element == "Al")
    return 1.26;

  if (element == "Si")
    return 1.16;

  if (element == "P")
    return 1.11;

  if (element == "S")
    return 1.03;

  if (element == "Cl")
    return 0.99;

  if (element == "Ar")
    return 0.96;

  /*4s*/
  if (element == "K")
    return 1.96;

  if (element == "Ca")
    return 1.71;

  /*3d*/
  if (element == "Sc")
    return 1.48;

  if (element == "Ti")
    return 1.36;

  if (element == "V")
    return 1.34;

  if (element == "Cr")
    return 1.22;

  if (element == "Mn")
    return 1.19;

  if (element == "Fe")
    return 1.16;

  if (element == "Co")
    return 1.11;

  if (element == "Ni")
    return 1.10;

  if (element == "Cu")
    return 1.12;

  if (element == "Zn")
    return 1.18;

  /*4p*/
  if (element == "Ga")
    return 1.24;

  if (element == "Ge")
    return 1.21;

  if (element == "As")
    return 1.21;

  if (element == "Se")
    return 1.16;

  if (element == "Br")
    return 1.14;

  if (element == "Kr")
    return 1.17;

  /*5s*/
  if (element == "Rb")
    return 2.10;

  if (element == "Sr")
    return 1.85;

  /*4d*/
  if (element == "Y")
    return 1.63;

  if (element == "Zr")
    return 1.54;

  if (element == "Nb")
    return 1.47;

  if (element == "Mo")
    return 1.38;

  if (element == "Tc")
    return 1.28;

  if (element == "Ru")
    return 1.25;

  if (element == "Rh")
    return 1.25;

  if (element == "Pd")
    return 1.20;

  if (element == "Ag")
    return 1.28;

  if (element == "Cd")
    return 1.36;

  /*5p*/
  if (element == "In")
    return 1.42;

  if (element == "Sn")
    return 1.40;

  if (element == "Sb")
    return 1.40;

  if (element == "Te")
    return 1.36;

  if (element == "I")
    return 1.33;

  if (element == "Xe")
    return 1.31;

  /*6s*/
  if (element == "Cs")
    return 2.32;

  if (element == "Ba")
    return 1.96;

  /*4f*/
  if (element == "Ce")
    return 1.63;

  if (element == "Pr")
    return 1.76;

  if (element == "Nd")
    return 1.74;

  if (element == "Pm")
    return 1.73;

  if (element == "Sm")
    return 1.72;

  if (element == "Eu")
    return 1.68;

  if (element == "Gd")
    return 1.69;

  if (element == "Tb")
    return 1.68;

  if (element == "Dy")
    return 1.67;

  if (element == "Ho")
    return 1.66;

  if (element == "Er")
    return 1.65;

  if (element == "Tm")
    return 1.64;

  if (element == "Yb")
    return 1.70;

  if (element == "Lu")
    return 1.62;

  /*5d*/
  if (element == "La")
    return 1.80;

  if (element == "Hf")
    return 1.52;

  if (element == "Ta")
    return 1.46;

  if (element == "W")
    return 1.37;

  if (element == "Re")
    return 1.31;

  if (element == "Os")
    return 1.29;

  if (element == "Ir")
    return 1.22;

  if (element == "Pt")
    return 1.23;

  if (element == "Au")
    return 1.24;

  if (element == "Hg")
    return 1.33;

  /*6p*/
  if (element == "Tl")
    return 1.44;

  if (element == "Pb")
    return 1.44;

  if (element == "Bi")
    return 1.51;

  if (element == "Po")
    return 1.45;

  if (element == "At")
    return 1.47;

  if (element == "Rn")
    return 1.42;

  /*7s*/
  if (element == "Fr")
    return 2.23;

  if (element == "Ra")
    return 2.01;

  /*5f*/
  if (element == "Th")
    return 1.75;

  if (element == "Pa")
    return 1.69;

  if (element == "U")
    return 1.70;

  if (element == "Np")
    return 1.71;

  if (element == "Pu")
    return 1.72;

  if (element == "Am")
    return 1.66;

  if (element == "Cm")
    return 1.66;

  if (element == "Bk")
    return 1.68;

  if (element == "Cf")
    return 1.68;

  if (element == "Es")
    return 1.65;

  if (element == "Fm")
    return 1.67;

  if (element == "Md")
    return 1.73;

  if (element == "No")
    return 1.76;

  if (element == "Lr")
    return 1.61;

  /*6d*/
  if (element == "Ac")
    return 1.86;

  if (element == "Rf")
    return 1.57;

  if (element == "Db")
    return 1.49;

  if (element == "Sg")
    return 1.43;

  if (element == "Bh")
    return 1.41;

  if (element == "Hs")
    return 1.34;

  if (element == "Mt")
    return 1.29;

  if (element == "Ds")
    return 1.28;

  if (element == "Rg")
    return 1.21;

  if (element == "Cn")
    return 1.22;

  /*7p*/
  if (element == "Nh")
    return 1.36;

  if (element == "Fl")
    return 1.43;

  if (element == "Mc")
    return 1.62;

  if (element == "Lv")
    return 1.75;

  if (element == "Ts")
    return 1.65;

  if (element == "Og")
    return 1.57;

  return 1.0; // Guard radius
} // covalentRadii

/*-----------------------------------------------------------------------------
 * Atomic form factors parameters are adimensional.
 * The atomic form factor parameters were taken from the International Tables
 * for Crystallography: http://it.iucr.org/CB/ch6o1v0001/
 *-----------------------------------------------------------------------------
 */
std::vector<double> atomicFormFactorParameters(std::string_view element) {
  /*1s*/
  if (element == "H")
    return {0.489918, 20.6593,  0.262003, 7.74039, 0.196767,
            49.5519,  0.049879, 2.20159,  0.001305};

  if (element == "He")
    return {0.8734,  9.1037, 0.6309, 3.3568, 0.3112,
            22.9276, 0.178,  0.9821, 0.0064};

  /*2s*/
  if (element == "Li")
    return {1.1282,  3.9546, 0.7508,  1.0524, 0.6175,
            85.3905, 0.4653, 168.261, 0.0377};

  if (element == "Be")
    return {1.5919,  43.6427, 1.1278, 1.8623, 0.5391,
            103.483, 0.7029,  0.542,  0.0385};

  /*2p*/
  if (element == "B")
    return {2.0545,  23.2185, 1.3326, 1.021,  1.0979,
            60.3498, 0.7068,  0.1403, -0.1932};

  if (element == "C")
    return {2.31,   20.8439, 1.02,    10.2075, 1.5886,
            0.5687, 0.865,   51.6512, 0.2156};

  if (element == "N")
    return {12.2126, 0.0057, 3.1322, 9.8933, 2.0125,
            28.9975, 1.1663, 0.5826, -11.529};

  if (element == "O")
    return {3.0485, 13.2771, 2.2868,  5.7011, 1.5463,
            0.3239, 0.867,   32.9089, 0.2508};

  if (element == "F")
    return {3.5392, 10.2825, 2.6412,  4.2944, 1.517,
            0.2615, 1.0243,  26.1476, 0.2776};

  if (element == "Ne")
    return {3.9553, 8.4042, 3.1125,  3.4262, 1.4546,
            0.2306, 1.1251, 21.7184, 0.3515};

  /*3s*/
  if (element == "Na")
    return {4.7626, 3.285,  3.1736,  8.8422, 1.2674,
            0.3136, 1.1128, 129.424, 0.676};

  if (element == "Mg")
    return {5.4204, 2.8275, 2.1735, 79.2611, 1.2269,
            0.3808, 2.3073, 7.1937, 0.8584};

  /*3p*/
  if (element == "Al")
    return {6.4202,  3.0387, 1.9002,  0.7426, 1.5936,
            31.5472, 1.9646, 85.0886, 1.1151};

  if (element == "Si")
    return {6.2915, 2.4386, 3.0353,  32.3337, 1.9891,
            0.6785, 1.541,  81.6937, 1.1407};

  if (element == "P")
    return {6.4345, 1.9067, 4.1791,  27.157, 1.78,
            0.526,  1.4908, 68.1645, 1.1149};

  if (element == "S")
    return {6.9053, 1.4679, 5.2034, 22.2151, 1.4379,
            0.2536, 1.5863, 56.172, 0.8669};

  if (element == "Cl")
    return {11.4604, 0.0104, 7.1964,  1.1662, 6.2556,
            18.5194, 1.6455, 47.7784, -9.5574};

  if (element == "Ar")
    return {7.4845,  0.9072, 6.7723,  14.8407, 0.6539,
            43.8983, 1.6442, 33.3929, 1.4445};

  /*4s*/
  if (element == "K")
    return {8.2186,  12.7949, 7.4398,  0.7748, 1.0519,
            213.187, 0.8659,  41.6841, 1.4228};

  if (element == "Ca")
    return {8.6266,  10.4421, 7.3873,  0.6599, 1.5899,
            85.7484, 1.0211,  178.437, 1.3751};

  /*3d*/
  if (element == "Sc")
    return {9.189,   9.0213, 7.3679,  0.5729, 1.6409,
            136.108, 1.468,  51.3531, 1.3329};

  if (element == "Ti")
    return {9.7595,  7.8508, 7.3558,  0.5,   1.6991,
            35.6338, 1.9021, 116.105, 1.2807};

  if (element == "V")
    return {10.2971, 6.8657, 7.3511,  0.4385, 2.0703,
            26.8938, 2.0571, 102.478, 1.2199};

  if (element == "Cr")
    return {10.6406, 6.1038, 7.3537,  0.392, 3.324,
            20.2626, 1.4922, 98.7399, 1.1832};

  if (element == "Mn")
    return {11.2819, 5.3409, 7.3573,  0.3432, 3.0193,
            17.8674, 2.2441, 83.7543, 1.0896};

  if (element == "Fe")
    return {11.7695, 4.7611, 7.3573,  0.3072, 3.5222,
            15.3535, 2.3045, 76.8805, 1.0369};

  if (element == "Co")
    return {12.2841, 4.2791, 7.3409,  0.2784, 4.0034,
            13.5359, 2.3488, 71.1692, 1.0118};

  if (element == "Ni")
    return {12.8376, 3.8785, 7.292,   0.2565, 4.4438,
            12.1763, 2.38,   66.3421, 1.0341};

  if (element == "Cu")
    return {13.338,  3.5828, 7.1676,  0.247, 5.6158,
            11.3966, 1.6735, 64.8126, 1.191};

  if (element == "Zn")
    return {14.0743, 3.2655, 7.0318,  0.2333, 5.1652,
            10.3163, 2.41,   58.7097, 1.3041};

  /*4p*/
  if (element == "Ga")
    return {15.2354, 3.0669, 6.7006,  0.2412, 4.3591,
            10.7805, 2.9623, 61.4135, 1.7189};

  if (element == "Ge")
    return {16.0816, 2.8509, 6.3747,  0.2516, 3.7068,
            11.4468, 3.683,  54.7625, 2.1313};

  if (element == "As")
    return {16.6723, 2.6345, 6.0701,  0.2647, 3.4313,
            12.9479, 4.2779, 47.7972, 2.531};

  if (element == "Se")
    return {17.0006, 2.4098, 5.8196,  0.2726, 3.9731,
            15.2372, 4.3543, 43.8163, 2.8409};

  if (element == "Br")
    return {17.1789, 2.1723, 5.2358,  16.5796, 5.6377,
            0.2609,  3.9851, 41.4328, 2.9557};

  if (element == "Kr")
    return {17.3555, 1.9384, 6.7286,  16.5623, 5.5493,
            0.2261,  3.5375, 39.3972, 2.825};

  /*5s*/
  if (element == "Rb")
    return {17.1784, 1.7888, 9.6435,  17.3151, 5.1399,
            0.2748,  1.5292, 164.934, 3.4873};

  if (element == "Sr")
    return {17.5663, 1.5564, 9.8184,  14.0988, 5.422,
            0.1664,  2.6694, 132.376, 2.5064};

  /*4d*/
  if (element == "Y")
    return {17.776,   1.4029,  10.2946, 12.8006, 5.72629,
            0.125599, 3.26588, 104.354, 1.91213};

  if (element == "Zr")
    return {17.8765,  1.27618, 10.948,  11.916, 5.41732,
            0.117622, 3.65721, 87.6627, 2.06929};

  if (element == "Nb")
    return {17.6142,  1.18865, 12.0144, 11.766, 4.04183,
            0.204785, 3.53346, 69.7957, 3.75591};

  if (element == "Mo")
    return {3.7025, 0.2772, 17.2356, 1.0958, 12.8876,
            11.004, 3.7429, 61.6584, 4.3875};

  if (element == "Tc")
    return {19.1301, 0.864132, 11.0948, 8.14487, 4.64901,
            21.5707, 2.71263,  86.8472, 5.40428};

  if (element == "Ru")
    return {19.2674, 0.80852, 12.9182, 8.43467, 4.86337,
            24.7997, 1.56756, 94.2928, 5.37874};

  if (element == "Rh")
    return {19.2957, 0.751536, 14.3501, 8.21758, 4.73425,
            25.8749, 1.28918,  98.6062, 5.328};

  if (element == "Pd")
    return {19.3319, 0.698655, 15.5017, 7.98929, 5.29537,
            25.2052, 0.605844, 76.8986, 5.26593};

  if (element == "Ag")
    return {19.2808, 0.6446, 16.6885, 7.4726, 4.8045,
            24.6605, 1.0463, 99.8156, 5.179};

  if (element == "Cd")
    return {19.2214, 0.5946, 17.6444, 6.9089, 4.461,
            24.7008, 1.6029, 87.4825, 5.0694};

  /*5p*/
  if (element == "In")
    return {19.1624, 0.5476, 18.5596, 6.3776, 4.2948,
            25.8499, 2.0396, 92.8029, 4.9391};

  if (element == "Sn")
    return {19.1889, 5.8303, 19.1005, 0.5031, 4.4585,
            26.8909, 2.4663, 83.9571, 4.7821};

  if (element == "Sb")
    return {19.6418, 5.3034, 19.0455, 0.4607, 5.0371,
            27.9074, 2.6827, 75.2825, 4.5909};

  if (element == "Te")
    return {19.9644, 4.81742, 19.0138, 0.420885, 6.14487,
            28.5284, 2.5239,  70.8403, 4.352};

  if (element == "I")
    return {20.1472, 4.347,  18.9949, 0.3814, 7.5138,
            27.766,  2.2735, 66.8776, 4.0712};

  if (element == "Xe")
    return {20.2933, 3.9282, 19.0298, 0.344, 8.9767,
            26.4659, 1.99,   64.2658, 3.7118};

  /*6s*/
  if (element == "Cs")
    return {20.3892, 3.569,  19.1062, 0.3107, 10.662,
            24.3879, 1.4953, 213.904, 3.3352};

  if (element == "Ba")
    return {20.3361, 3.216,  19.297,  0.2756, 10.888,
            20.2073, 2.6959, 167.202, 2.7731};

  /*4f*/
  if (element == "Ce")
    return {21.1671, 2.81219, 19.7695, 0.226836, 11.8513,
            17.6083, 3.33049, 127.113, 1.86264};

  if (element == "Pr")
    return {22.044,  2.77393, 19.6697, 0.222087, 12.3856,
            16.7669, 2.82428, 143.644, 2.0583};

  if (element == "Nd")
    return {22.6845, 2.66248, 19.6847, 0.210628, 12.774,
            15.885,  2.85137, 137.903, 1.98486};

  if (element == "Pm")
    return {23.3405, 2.5627,  19.6095, 0.202088, 13.1235,
            15.1009, 2.87516, 132.721, 2.02876};

  if (element == "Sm")
    return {24.0042, 2.47274, 19.4258, 0.196451, 13.4396,
            14.3996, 2.89604, 128.007, 2.20963};

  if (element == "Eu")
    return {24.6274, 2.3879, 19.0886, 0.1942, 13.7603,
            13.7546, 2.9227, 123.174, 2.5745};

  if (element == "Gd")
    return {25.0709, 2.25341, 19.0798, 0.181951, 13.8518,
            12.9331, 3.54545, 101.398, 2.4196};

  if (element == "Tb")
    return {25.8976, 2.24256, 18.2185, 0.196143, 14.3167,
            12.6648, 2.95354, 115.362, 3.58324};

  if (element == "Dy")
    return {26.507,  2.1802,  17.6383, 0.202172, 14.5596,
            12.1899, 2.96577, 111.874, 4.29728};

  if (element == "Ho")
    return {26.9049, 2.07051, 17.294,  0.19794, 14.5583,
            11.4407, 3.63837, 92.6566, 4.56796};

  if (element == "Er")
    return {27.6563, 2.07356, 16.4285, 0.223545, 14.9779,
            11.3604, 2.98233, 105.703, 5.92046};

  if (element == "Tm")
    return {28.1819, 2.02859, 15.8851, 0.238849, 15.1542,
            10.9975, 2.98706, 102.961, 6.75621};

  if (element == "Yb")
    return {28.6641, 1.9889,  15.4345, 0.257119, 15.3087,
            10.6647, 2.98963, 100.417, 7.56672};

  if (element == "Lu")
    return {28.9476,  1.90182, 15.2208, 9.98519, 15.1,
            0.261033, 3.71601, 84.3298, 7.97628};

  /*5d*/
  if (element == "La")
    return {20.578,  2.94817, 19.599,  0.244475, 11.3727,
            18.7726, 3.28719, 133.124, 2.14678};

  if (element == "Hf")
    return {29.144,   1.83262, 15.1726, 9.5999, 14.7586,
            0.275116, 4.30013, 72.029,  8.58154};

  if (element == "Ta")
    return {29.2024,  1.77333, 15.2293, 9.37046, 14.5135,
            0.295977, 4.76492, 63.3644, 9.24354};

  if (element == "W")
    return {29.0818,  1.72029, 15.43,  9.2259, 14.4327,
            0.321703, 5.11982, 57.056, 9.8875};

  if (element == "Re")
    return {28.7621, 1.67191, 15.7189, 9.09227, 14.5564,
            0.3505,  5.44174, 52.0861, 10.472};

  if (element == "Os")
    return {28.1894,  1.62903, 16.155,  8.97948, 14.9305,
            0.382661, 5.67589, 48.1647, 11.0005};

  if (element == "Ir")
    return {27.3049,  1.59279, 16.7296, 8.86553, 15.6115,
            0.417916, 5.83377, 45.0011, 11.4722};

  if (element == "Pt")
    return {27.0059,  1.51293, 17.7639, 8.81174, 15.7131,
            0.424593, 5.7837,  38.6103, 11.6883};

  if (element == "Au")
    return {16.8819, 0.4611, 18.5913, 8.6216, 25.5582,
            1.4826,  5.86,   36.3956, 12.0658};

  if (element == "Hg")
    return {20.6809, 0.545,  19.0417, 8.4484, 21.6575,
            1.5729,  5.9676, 38.3246, 12.6089};

  /*6p*/
  if (element == "Tl")
    return {27.5446, 0.65515, 19.1584, 8.70751, 15.538,
            1.96347, 5.52593, 45.8149, 13.1746};

  if (element == "Pb")
    return {31.0617, 0.6902, 13.0637, 2.3576, 18.442,
            8.618,   5.9696, 47.2579, 13.4118};

  if (element == "Bi")
    return {33.3689, 0.704,  12.951,  2.9238, 16.5877,
            8.7937,  6.4692, 48.0093, 13.5782};

  if (element == "Po")
    return {34.6726, 0.700999, 15.4733, 3.55078, 13.1138,
            9.55642, 7.02588,  47.0045, 13.677};

  if (element == "At")
    return {35.3163, 0.68587, 19.0211, 3.97458, 9.49887,
            11.3824, 7.42518, 45.4715, 13.7108};

  if (element == "Rn")
    return {35.5631, 0.6631, 21.2816, 4.0691, 8.0037,
            14.0422, 7.4433, 44.2473, 13.6905};

  /*7s*/
  if (element == "Fr")
    return {35.9299, 0.646453, 23.0547, 4.17619, 12.1439,
            23.1052, 2.11253,  150.645, 13.7247};

  if (element == "Ra")
    return {35.763,  0.616341, 22.9064, 3.87135, 12.4739,
            19.9887, 3.21097,  142.325, 13.6211};

  /*5f*/
  if (element == "Th")
    return {35.5645, 0.563359, 23.4219, 3.46204, 12.7473,
            17.8309, 4.80703,  99.1722, 13.4314};

  if (element == "Pa")
    return {35.8847, 0.547751, 23.2948, 3.41519, 14.1891,
            16.9235, 4.17287,  105.251, 13.4287};

  if (element == "U")
    return {36.0228, 0.5293, 23.4128, 3.3253, 14.9491,
            16.0927, 4.188,  100.613, 13.3966};

  if (element == "Np")
    return {36.1874, 0.511929, 23.5964, 3.25396, 15.6402,
            15.3622, 4.1855,   97.4908, 13.3573};

  if (element == "Pu")
    return {36.5254, 0.499384, 23.8083, 3.26371, 16.7707,
            14.9455, 3.47947,  105.98,  13.3812};

  if (element == "Am")
    return {36.6706, 0.483629, 24.0992, 3.20647, 17.3415,
            14.3136, 3.49331,  102.273, 13.3592};

  if (element == "Cm")
    return {36.6488, 0.465154, 24.4096, 3.08997, 17.399,
            13.4346, 4.21665,  88.4834, 13.2887};

  if (element == "Bk")
    return {36.7881, 0.451018, 24.7736, 3.04619, 17.8919,
            12.8946, 4.23284,  86.003,  13.2754};

  if (element == "Cf")
    return {36.9185, 0.437533, 25.1995, 3.00775, 18.3317,
            12.4044, 4.24391,  83.7881, 13.2674};

  /*6d*/
  if (element == "Ac")
    return {35.6597, 0.589092, 23.1032, 3.65155, 12.5977,
            18.599,  4.08655,  117.02,  13.5266};

  return {36.9185, 0.437533, 25.1995, 3.00775, 18.3317,
          12.4044, 4.24391,  83.7881, 13.2674}; // Guard parameters
}; // atomicFormFactorParameters
