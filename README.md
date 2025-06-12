# InitialIdealsRegularSubdivisions

This is the online supplement to the paper "Regular subdivisions, bounds on initial ideals, and categorical limits" (<a href="https://arxiv.org/abs/2411.12819">arXiv 2411.12819</a>) by George Balla, Daniel Corey, Igor Makhlin, and Victoria Schleis.

We use OSCAR version 1.3.1 in Julia. It may be essential to use this version of OSCAR. To ensure this, do the following. First, run `julia --project=.` in the terminal from the root of this project. Next, open `julia` and run the following: 

```
julia> using Pkg
julia> Pkg.instantiate()
```
All of our code is contained in the notebook `Computing_I_w_and_Omega.ipynb`. In particular, this notebook contains the details from the following parts of the paper.

**Example 5.4.** We compute $\Omega(I)$ and show that it is a union of 72 of the 102 maximal cones of $\mathrm{Sec} \mathcal{A}(I)$. 

**Example 5.5.** We compute the fan $\Omega(I)$ for the Pluecker ideal $I = I_{3,6}$ and show that it has 30 maximal cones of dimension 10 and 10 more maximal cones of dimension 8. This uses precomputed data stored in the file `TGr36.mrdi` which was provided to us by Michael Joswig. 

**Appendix B.** We create functions that compute the point configuration $\mathcal{A}(I)$ associated to the ideal $I$. 

