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
  - name: Ulises Santiago
    affiliation: 3
  - name: Ariel A. Valladares
    affiliation: 1
affiliations:
 - name: Instituto de Investigaciones en Materiales, Universidad Nacional Autónoma de México
   index: 1
 - name: Facultad de Ciencias, Universidad Nacional Autónoma de México
   index: 2
 - name: Department of Physics and Astronomy, The University of Texas at San Antonio
   index: 3
date: 30 November 2020
bibliography: paper.bib
---

# Summary

For almost a century since J.D. Bernal attempt at molecular theory of liquid structure [@bernal_attempt_1937], correlation functions have been the bridge to compare theoretical calculations with experimental measurements in the study of disordered materials.

Pair-Distribution Function($g(r)$), Radial Distribution Function ($J(r)$), Plane-Angle Distribution ($g(\theta)$) and Coordination Numbers ($n_c$) have been widely used to characterize amorphous and liquid materials [@waseda_structure_1980; @elliott_physics_1986; @valladares_new_2011] and, in particular, Bulk Metallic Glasses [@miller_bulk_2007; @galvan-colin_short-range_2015].

**Correlation** is an Open-Source software designed to analyze liquid structures and amorphous solids; the software is user-friendly, the modular design make it easy to integrate in High-throughput computing (HTC) to process structures with a large number of constituents in a standardized fashion. **Correlation** is ready to be used in Windows, Linux and Mac. Currently, we support DMol3 (CAR), CASTEP and ONETEP (CELL), LAMMPS (XYZ) and VASP (OUTCAR) structure files. The code can handle up to 25,000 atoms, so it can be used to analyze both classical simulations and first-principles simulations. At the end, the output of every single correlation function is exported in it's own individual comma-separated value file (CSV), to further analyze the results.

# Statement of Need

As time goes by the number of atoms in theoretical calculations has grown from a few dozens to hundreds of thousands of atoms, and with this increment the complexity to calculate the correlation functions that represents the structure of materials has steadily increase.
To answer this need, there has been several tools developed to calculate some of the most used correlations functions like: pair correlation function, radial distribution function and plane angle distributions.
However, the use of these tools has been limited, either by a prohibiting cost, been restricted to private academic groups, or geopolitical limitations introduced by the licensing, or being specialized to an specific material simulation software.

With these limitations in mind, we decided to create a software that could calculate the correlation functions of materials, as well as the more interesting properties derived from these functions. While making the software accessible to as many people as we could.

# Mathematical Background

## Pair Correlation Functions

The structure factor is one of the most useful tools in analyzing the scattering patterns obtained from X-ray, electron and neutron diffraction experiments of crystalline, disordered and amorphous materials.
The Fourier transform of the scattering intensity given by the structure factor $S(Q)$, yields the pair-distribution function $g(r)$ defined by:

\begin{equation}\label{eq:PDF}
G(r) = 4\pi\rho_0(g(r)-1)=\frac{2}{\pi}\int_{0}^{\inf} Q[S(Q)-1]sin(Qr)dQ,
\end{equation}
where $G(r)$ is the reduced pair-distribution function.

![Schematic depiction of the first and second neighbors coordination spheres for an amorphous metallic alloy and the corresponding pair-distribution function. \label{fig:RDF}](./Images/Fig1.png){ width=75% }

The pair-distribution function (PDF) could also be seen like a distance map inside the material, the $g(r)$ function gives the probability of finding two atoms separated by the distance ($r$) as can be seen in \autoref{fig:RDF}.

The PDF is normalized so that, as $r \to \infty$, $g(r) \to 1$, Also, as $r \to 0$ (for $r$ shorter than the distance of the closest approach of pair of atoms), $g(r)$ become zero. The main advantage of the PDF, and the related functions, is that they emphasize the short-range order of a material.

### Reduced pair-distribution function $G(r)$

One of the most widely used pair correlation function is the reduced pair-distribution function. This is defined as $G(r) = 4\pi \rho_0 (g(r)-1)$. From this definition, and the previously discussed tendency at large $r$ of the PDF, it's clear that the reduced pair-distribution function (rPDF) oscillates around zero as $r \to \infty$. It also becomes evident that as $r \to 0$ (for $r$ smaller than the closest pair of atoms), the rPDF behaves like $-4\pi \rho_0$.

While the PDF ($g(r)$) has an intuitive geometric definition, the rPDF ($G(r)$) can be directly obtained by a Fourier transformation of the structure factor ($S(Q)$) as can be seen in \autoref(eq:PDF), this close relation with the experimental data explains the popularity that this function has. Also, it has an other advantage over other correlation functions, like PDF ($g(r)$) as the numerical density $\rho_0 = N/V$ should be estimated to normalize the other correlation functions. This is not necessary in rPDF ($G(r)$), on the contrary this information is already contained in $G(r)$ as the slope of the function as $r \to 0$.

All of these advantages for the analysis of experimental measurements, are not true for theoretical simulations, as the numerical density ($\rho_0$) could be directly obtained by a simple division of the number of the particles and the volume of the cell, while the mapping of the $Q$ vector in the structure factor ($S(Q)$) increase in difficulty as the cell volume increase, as it's the case in disordered, amorphous and liquids.

### Radial distribution function $J(r)$

The last pair correlation function we discuss is also one of the most physical intuitive, The PDF, $g(r)$ is related to de Radial Distribution Function (RDF) by:

\begin{equation}\label{eq:RDF}
J(r) = 4\pi r^{2} \rho_0 g(r).
\end{equation}

The RDF has the useful property that the quantity $J(r)dr$ gives the number of atoms of an spherical shell with inner radius $r$ and a thickness of $dr$ around every atom as depicted in \autoref(fig:RDF). For example the coordination number, or the number of neighbors ($CN$), is given by:

\begin{equation}\label{eq:RDF}
CN = \int_{r_1}^{r_2} J(r) dr,
\end{equation}

where $r_1$ and $r_2$ define the RDF peak corresponding to the coordination shell in question.  

## Plane-Angle Distribution

The use of higher order correlation functions to analyze the structure of liquids and amorphous solids has been proposed in the literature [@hafner_triplet_1982, @galvan-colin_ab_2016], trying to reproduce the success obtained by J.D. Bernal in the analysis of the structure of liquids [@bernal_bakerian_1964].

In particular, the Bond-Angle Distribution $f(\theta)$ has been used to characterize the short-range order of amorphous and liquid structures[@galvan-colin_short-range_2015, @galvan-colin_ab_2016, @mata-pinzon_superconductivity_2016]. We propose to substitute the term “Bond-Angle Distribution” by the also frequently used “Plane-Angle Distribution” since in condensed matter, proximity does not necessarily imply bonding.

# Benchmarks

![Pair Distribution Functions $g(r)$ on the left, Plane-Angle Distributions on the right for: crystalline silicon, Graphene Layer, amorphous Palladim, amourphous Palladium Hydride and liquid Bismuth. Correlation in gray, Forcite in black, Plane-Angles in red \label{fig:RDF-PAD}](./Images/Fig2.png){ width=70% }


In order to assess the performance of **Correlation**, we calculate the PDF and PAD for two well known structures (Crystalline Silicon and a Graphene Layer), and compare the results with the commercially available software Forcite included in the Materials Studio suite [@materials_2016], to test **Correlation** in Amorphous and Liquids materials we selected amorphous palladim [@rodriguez_emergence_2019], amorphous palladium hydride [rodriguez_calculo_2019] and liquid Bismuth. Because of the complexity to calculate PAD of amorphous and liquids in Forcite, we choose to compare them with the code developed by U. Santiago as part of his PhD thesis in the workgroup of Prof. Valladares [@santiago_simulacion_2011].


The results of these benchmarks are shown in Figure \autoref{fig:RDF-PAD}, and the structures used to calculate these figures are included in the code as test 1 to 5. The last estructure included as test 6 is a 2x2x2 supercell of amorphous palladium hydride included in test 4, to benchmark Memory and CPU performance in an structure with thousands of atoms.

# Conclusion & Perspective

**Correlation** is a lightweight, modular software that can be use in HPC and adapted to analyze the main correlation functions used to characterize: crystalline, amorphous and liquids.

**Correlation** software has been used in previously published work [@galvan-colin_short-range_2015, @galvan-colin_ab_2016, @mata-pinzon_superconductivity_2016] as well as several PhD Thesis [@santiago_simulacion_2011, @romero-rangel_simulaciones_2014, @mejia-mendoza_estudio_2014, @galvan-colin_atomic_2016, @mata-pinzon_propiedades_2016, @rodriguez_calculo_2019]. And our workgroup will continue to support and enrich the software for the foreseeable future. We are open to receive suggestions that further improve the functionality of the software.

# Acknowledgements

I.R. acknowledge PAPIIT, DGAPA-UNAM for his posdoctoral fellowship.
D.H.R. acknowledge Consejo Nacional de Ciencia y Tecnología (CONACyT) for supporting his graduate studies.
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under Grant No. IN104617 and IN116520.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out in the Supercomputing Center of DGTIC-UNAM.
I.R. would like to express his gratitude to F. B. Quiroga, M. A. Carrillo, R. S. Vilchis, S. Villareal and A. de Leon, for their time invested in testing the code, as well as the structures provided for benchmarks and tests.

# References
