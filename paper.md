---
title: 'Correlation: An Analyzing Tool for Liquids and for Amorphous Solids'
tags:
  - C++
  - Materials
  - Molecular dynamics
  - Pair Correlation Function
  - Plane Angle Distribution
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
 - name: Department of Computational and Systems Biology, University of Pittsburgh
   index: 3
date: 30 November 2020
bibliography: paper.bib
csl: aps.csl
---

# Summary

For almost a century, since Bernal's attempts at a molecular theory of liquid structure [@bernal_attempt_1937], correlation functions have been the bridge to compare theoretical calculations with experimental measurements in the study of disordered materials.

Pair Distribution Functions ($g(r)$), Radial Distribution Functions ($J(r)$), Plane Angle Distributions ($g(\theta)$) and Coordination Numbers ($n_c$) have been widely used to characterize amorphous and liquid materials [@waseda_structure_1980; @elliott_physics_1986; @valladares_new_2011] and, in particular Bulk Metallic Glasses [@miller_bulk_2007; @galvan-colin_short-range_2015].

**Correlation** is an Open-Source software designed to analyze liquid structures and amorphous solids; the software is user-friendly, the modular design makes it easy to integrate in High-Throughput Computing (HTC) to process structures with a large number of constituents in a standardized fashion. **Correlation** is ready to be used in Windows, Linux and Mac. Currently, we support DMol3 (CAR), CASTEP (CELL), ONETEP (DAT) and VASP (POSCAR) structure files. The code can handle up to 100,000 atoms, so it can be used to analyze both classical and first-principles simulations. At the end, the output of every single correlation function is exported to the corresponding comma-separated value file (CSV), to further analyze the results.

# Statement of Need

As time goes by, the number of atoms in theoretical calculations has grown from a few dozens to hundreds of thousands of atoms, and with this increment the complexity to calculate the correlation functions that represent the structure of materials has steadily increased.
To answer this need, there have been several tools developed to calculate some of the most used correlation functions like: pair distribution functions, radial distribution functions, and plane angle distributions, here we present an incomplete list of these tools:
-	**Forcite Plus**: Forcite Plus is part of the Materials Studio suite [@materials_2016], the program includes analyzing tools to compute structure properties like RDF and PAD.
-	**PTRAJ/CPPTRAJ**: A tool designed to analyze Amber Molecular Dynamics trajectories and related properties including RDF and time correlation functions, included in the AmbarTools suite [@roe_ptraj_2013].
-	**VASPKIT**: A post-processing tool for VASP calculated data [@VASPKIT-2021], the code includes tools to analyze structural properties and dynamics trajectories.
-	[**rdfpy**](https://github.com/by256/rdfpy): An open Python library for fast computation of 2D and 3D radial distribution functions of crystalline structures in the Crystallographic Information File (CIF).
-	[**RadialDistributionFunction **](https://github.com/RaulPPelaez/RadialDistributionFunction): Computes the Radial Distribution Function (RDF) of a group of positions in a file, averages it for all snapshots in the file, atom positions must be in a custom format.

However, the use of these tools has been limited, either by a prohibiting cost (**Forcite Plus**), or has been restricted to private academic groups, or geopolitical limitations introduced by the licensing (**CASTEP** postprocessing tools), or by being specially designed to specific software for material simulation (**PTRAJ/CPPTRAJ**, **VASPKIT**), or by having a narrow scope of input formats and correlation functions calculated (**rdfpy**, **RadialDistributionFunction **).
With these limitations in mind, we decided to create software that could calculate several correlation functions of materials, as well as the more interesting properties derived from these functions. While making the software accessible to as many people as possible.

# Mathematical Background

## Pair Distribution Functions

The structure factor is one of the most useful tools in analyzing the scattering patterns obtained from X-ray, electron and neutron diffraction experiments of crystalline, disordered and amorphous materials.
The Fourier transform of the scattering intensity given by the structure factor $S(Q)$, yields the pair distribution function (PDF) $g(r)$ defined by:

\begin{equation}\label{eq:PDF}
G(r) = 4\pi\rho_0(g(r)-1)=\frac{2}{\pi}\int_{0}^{\infty} Q[S(Q)-1]sin(Qr)dQ,
\end{equation}
where $G(r)$ is the reduced pair distribution function (rPDF) which will be discussed next.

![Schematic depiction of the first and second neighbor’s coordination spheres for an amorphous metallic alloy and the corresponding pair distribution function. Design inspired by J. M. Ziman, Models of disorder (Cambridge University Press, 1979). \label{fig:RDF}](./Images/Fig1.png){ width=75% }

The pair distribution function could also be seen like a distance map inside the material, the $g(r)$ function gives how feasible is finding two atoms separated by the distance ($r$) as can be seen in \autoref{fig:RDF} [@ziman_models_1979].

The PDF is normalized so that, as $r \to \infty$, $g(r) \to 1$, Also, as $r \to 0$ (for $r$ shorter than the distance of the closest approach of pair of atoms), $g(r)$ becomes zero. The main advantage of the PDF and the related functions, is that they emphasize the short-range order of a material.

### Reduced pair distribution function $G(r)$

One of the most widely used pair correlation function is the reduced pair distribution function. This is defined as $G(r) = 4\pi \rho_0 (g(r)-1)$. From this definition, and the previously discussed tendency at large $r$ of the PDF, it's clear that the reduced pair distribution function (rPDF) oscillates around zero as $r \to \infty$. It also becomes evident that as $r \to 0$ (for $r$ smaller than the closest pair of atoms), the rPDF behaves like $-4\pi \rho_0$.

While the PDF ($g(r)$) has an intuitive geometric definition, the rPDF ($G(r)$) can be directly obtained by a Fourier transform of the structure factor ($S(Q)$) as can be seen in \autoref{eq:PDF}; this close relation with the experimental data explains the popularity that this function has. Also, it has another advantage over other correlation functions [like PDF ($g(r)$)] as the numerical density $\rho_0 = N/V$ needs to be estimated in order to normalize the functions. This is not necessary in rPDF ($G(r)$); where this information is already contained in $G(r)$ as the slope of the function when $r \to 0$.

### Radial distribution function $J(r)$

The last correlation function we shall discuss is also one of the most physically intuitive, the RDF, $J(r)$, which is related to the Pair Distribution Function (PDF) by:

\begin{equation}\label{eq:RDF}
J(r) = 4\pi r^{2} \rho_0 g(r).
\end{equation}

The RDF has the useful property that the quantity $J(r)dr$ gives the number of atoms in a spherical shell with inner radius $r$ and thickness $dr$ around every atom as depicted in \autoref{fig:RDF}. For example, the coordination number, or the number of neighbors ($n_c$), is given by:

\begin{equation}\label{eq:CN}
n_c = \int_{r_1}^{r_2} J(r) dr,
\end{equation}

where $r_1$ and $r_2$ define the RDF peak corresponding to the coordination shells in question.  

## Plane Angle Distribution

The use of higher order correlation functions to analyze the structure of liquids and amorphous solids has been proposed in the literature [@hafner_triplet_1982; @galvan-colin_ab_2016], trying to reproduce the success obtained by Bernal in the analysis of the structure of liquids [@bernal_bakerian_1964].

In particular, the Plane Angle Distribution, also known as the Bond Angle Distribution $f(\theta)$ has been used to characterize the short-range order of amorphous and liquid structures [@galvan-colin_short-range_2015; @galvan-colin_ab_2016; @mata-pinzon_superconductivity_2016]. Here we propose to substitute the term “Bond Angle Distribution” (BAD), by the term “Plane Angle Distribution” (PAD), also frequently used, since in condensed matter, proximity does not necessarily imply bonding.

# Benchmarks

![Pair Distribution Functions $g(r)$ on the left, Plane Angle Distributions on the right for: crystalline silicon, graphene Layer, amorphous palladium, amorphous palladium hydride and liquid bismuth. Correlation in gray, Forcite in black, Plane Angles in red. Similarity is remarkable between **Correlation** and Forcite as can be seen in all PDFs. Figures (a) to (d) indicate that the coincidence in the two results overlap completely. \label{fig:RDF-PAD}](./Images/Fig2.png){ width=75% }


In order to assess the performance of **Correlation**, we calculated the PDF and PAD for two well known structures (Crystalline Silicon and a Graphene Layer), and compared the results with the commercially available software Forcite included in the Materials Studio suite [@materials_2016]; to test **Correlation** in amorphous and liquid materials we selected amorphous palladium [@rodriguez_emergence_2019], amorphous palladium hydride [@rodriguez_calculo_2019] and liquid bismuth. Because of the complexity to calculate PAD of amorphous and liquids in Forcite, we chose to compare them with the code developed by U. Santiago within our group. [@santiago_simulacion_2011].


The results of these benchmarks are shown in \autoref{fig:RDF-PAD}, and the structures used to calculate these figures are included in the code as tests 1 to 5. The last structure included as test 6 is a 2x2x2 supercell of amorphous palladium hydride included in test 4, to benchmark memory and CPU performance in a structure with thousands of atoms.

# Conclusion & Perspective

**Correlation** is a lightweight, modular software that can be used in HTC and adapted to analyze the main correlation functions used to characterize: crystalline, amorphous solids, and liquids.

**Correlation** software has been used in previously published work [@galvan-colin_short-range_2015; @galvan-colin_ab_2016; @mata-pinzon_superconductivity_2016] as well as several Ph.D. Theses [@santiago_simulacion_2011; @romero-rangel_simulaciones_2014; @mejia-mendoza_estudio_2014; @galvan-colin_atomic_2016; @mata-pinzon_propiedades_2016; @rodriguez_calculo_2019] developed in our group.
We will continue to support and enrich the software in the foreseeable future, here we list the features planned to be added in the future:
-	Support for different kinds of output files like hdf5 standard.
-	Inclusion of other correlation functions like Velocity Correlation Functions, to further improve the analysis of liquids and phase transitions.
-	Inclusion of structure factor and x-ray diffraction, to facilitate the comparison with experimental results.
-	Parallelization of the main loop, to further improve the code by switching to a ‘divide-and-conquer paradigm.

We are open to receive suggestions that would further improve the functionality of the software. Address all comments and observations to the GitHub: [Correlation](https://github.com/Isurwars/Correlation), or to the email of the first author, I.R: [isurwars@ciencias.unam.mx]( isurwars@ciencias.unam.mx)

# Acknowledgements

I.R. acknowledges PAPIIT, DGAPA-UNAM for his posdoctoral fellowship.
D.H.R. acknowledges Consejo Nacional de Ciencia y Tecnología (CONACyT) for supporting his graduate studies.
A.A.V., R.M.V., and A.V. thank DGAPA-UNAM for continued financial support to carry out research projects under grants No. IN104617 and IN116520.
M. T. Vázquez and O. Jiménez provided the information requested.
A. López and A. Pompa helped with the maintenance and support of the supercomputer in IIM-UNAM.
Simulations were partially carried out at the Supercomputing Center of DGTIC-UNAM.
We would like to express our gratitude to F. B. Quiroga, M. A. Carrillo, R. S. Vilchis, S. Villarreal and A. de León, for the time invested in testing the code, as well as the structures provided for benchmarks and tests.

# References
