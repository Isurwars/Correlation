title: ' Correlation: An Analyzing Tool for Liquids and for Amorphous Solids'
tags:
  - C++
  - Materials
  - Molecular dynamics
  - Pair Correlation Function
  - Plane-Angle Distribution
authors:
  - name: Isaías Rodríguez]
    orcid: 0000-0002-4359-2742
    affiliation: "1" # (Multiple affiliations must be quoted)
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


# Summary

Pair Correlation Functions (g(r)), Radial Distribution Functions (J(r)), Plane-Angle Distributions (g(θ)) and Coordination Numbers (CN) have been widely used to characterize liquid and amorphous materials [1-3] and, in particular, Bulk Metallic Glasses [4,5].

Correlation is an Open-Source software designed by our group to analyze liquid structures and amorphous solids; the software is user-friendly, and easy to integrate in High-throughput computing (HTC) to process structures with a large number of constituents in a standard fashion. We propose to substitute the term “Bond-Angle Distribution” by the also frequently used “Plane-Angle Distribution” since in condensed matter, proximity does not necessarily imply bonding. The software is ready to be used in Windows, Linux and Mac. Currently, we support DMol3 (CAR), CASTEP and ONETEP (CELL), LAMMPS (XYZ) and VASP (OUTCAR) structure files. The code can handle up to 25,000 atoms, so it can be used to analyze both classical simulations and first-principles simulations. At the end, the output is exported in comma-separated value files (CSV), to further analyze the results.

# Introduction

## Pair Correlation Function



Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
