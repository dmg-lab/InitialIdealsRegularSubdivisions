# InitialIdealsRegularSubdivisions

This repository contains the online supplement to the preprint "Regular subdivisions, bounds on initial ideals, and categorical limits" (<a href="https://arxiv.org/abs/2411.12819">arXiv 2411.12819</a>) by George Balla, Daniel Corey, Igor Makhlin, and Victoria Schleis, as well as several additional computations not currently discussed in the paper.

We use the Julia package OSCAR (version 1.6.0), the Polymake.jl package, and Gfan.

The notebook `Computing_I_w_and_Omega.ipynb` contains an overview of our main functions, which are also found in `InitialIdealsRegularSubdivisions.jl`. Furthermore,
 - the notebook `omega_Gr2n.ipynb` contains computations of the full set $\Omega(I_{2,n})$ for small $n$;
 - the notebook `schubert.ipynb` computes the restricted sets $\Omega_t(I)$ for Schubert varieties in Grassmannians;
 - the notebook `spinor.ipynb` deals with spinor varieties.

The file `gfan_functions.jl` contains an interface for converting data to and from Gfan. The directory contains precomputations of large tropicalizations and secondary fans.

