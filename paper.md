---
title: ' Correlation: An Analyzing Tool for Liquids and for Amorphous Solids'
tags:
  - C++
  - Materials
  - Molecular dynamics
  - Pair Correlation Function
  - Plane-Angle Distribution
authors:
  - name: Isaías Rodríguez
    orcid: 0000-0002-4359-2742
    affiliation: 1
  - name: Renela M. Valladares
    affiliation: 2
  - name: Alexander Valladares
    affiliation: 2
  - name: David Hinojosa-Romero
    affiliation: 1
  - name: Ariel A. Valladares
    affiliation: 1
affiliations:
 - name: Instituto de Investigaciones en Materiales, Universidad Nacional Autónoma de México
   index: 1
 - name: Facultad de Ciencias, Universidad Nacional Autónoma de México
   index: 2
date: 30 November 2020
bibliography: paper.bib
---

# Summary

Pair Correlation Functions ($g(r)$), Radial Distribution Functions ($J(r)$), Plane-Angle Distributions ($g(\theta)$) and Coordination Numbers (CN) have been widely used to characterize liquid and amorphous materials [@waseda_structure_1980; @elliott_physics_1986; @valladares_new_2011] and, in particular, Bulk Metallic Glasses [@miller_bulk_2007; @galvan-colin_short-range_2015].

Correlation is an Open-Source software designed by our group to analyze liquid structures and amorphous solids; the software is user-friendly, and easy to integrate in High-throughput computing (HTC) to process structures with a large number of constituents in a standard fashion. We propose to substitute the term “Bond-Angle Distribution” by the also frequently used “Plane-Angle Distribution” since in condensed matter, proximity does not necessarily imply bonding. The software is ready to be used in Windows, Linux and Mac. Currently, we support DMol3 (CAR), CASTEP and ONETEP (CELL), LAMMPS (XYZ) and VASP (OUTCAR) structure files. The code can handle up to 25,000 atoms, so it can be used to analyze both classical simulations and first-principles simulations. At the end, the output is exported in comma-separated value files (CSV), to further analyze the results.

# Introduction

For almost a century since J.D. Bernal attempt at molecular theory of liquid structure [@bernal_attempt_1937], correlation functions have been the bridge to compare theoretical calculations with experimental measurements in the study of disordered materials.

As time goes by the number of atoms in theoretical calculations has grown from a few dozens to hundreds of thousands of atoms, and with this increment the complexity to calculate the correlation functions that represents the structure of materials has steadily increase.
To answer this need, there has been several tools developed to calculate some of the most used correlations functions like: pair correlation function, radial distribution function and plane angle distributions.
However, the use of these tools has been limited, either by a prohibiting cost, been restricted to private academic groups, or geopolitical limitations introduced by the licensing, or being specific to an specific software.

With these limitations in mind, we decided to create a software that could calculate the correlation functions of materials, as well as the more interesting properties derived from these functions. While making the software accessible to as many people as we could.


# Theory

## Pair Correlation Functions

The structure factor is one of the most useful tools in analyzing the scattering patterns obtained from X-ray, electron and neutron diffraction experiments of crystalline, disordered and amorphous materials.
The Fourier transform of the scattering intensity given by the structure factor $S(Q)$, yields the pair-distribution function $g(r)$ defined by:
$$
G(r) = 4\pi\rho_0(g(r)-1)=\frac{2}{\pi}\int_{0}^{\inf} Q[S(Q)-1]sin(Qr)dQ,
$$
where $G(r)$ is the reduced pair-distribution function.

![Schematic PDF \label{fig:RDF}](./Images/Fig1.png){ width=75% }

**Figure 1:** Schematic depiction of the first and second neighbors coordination spheres for an amorphous metallic alloy and the corresponding pair-distribution function $g(r)$.

The pair-distribution function (PDF) could also be seen like a distance map inside the material, the $g(r)$ function gives the probability of finding two atoms separated by the distance ($r$) as can be seen in \autoref{fig:RDF}.


# Acknowledgements

I.R. acknowledge PAPIIT, DGAPA-UNAM for his posdoctoral fellowship.
D.H.R. acknowledge Consejo Nacional de Ciencia y Tecnología (CONACyT) for supporting his graduate studies.
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under Grant No. IN104617 and IN116520.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out in the Supercomputing Center of DGTIC-UNAM.

# References
