# Toyfold 3
Calculations and simulations of RNA tertiary structure in 3D
_Author: Rhiju Das, 2020_

## Goals
* Allow rapid simulations to test partitition function calcs -- including tertiary structure -- in MATLAB.
* Actual predictions for `C_eff` and structural ensembles for circular RNA
* Actual predictions for nearest neighbor energies and structural ensembles for helices
* Test factorization of free energy in terms of local motif energy (including K_d); motif modularity costs; tertiary closure.
* Test simple analytic representations of translation/rotation SE(3) distributions in tertiary closure costs, including Gaussian models and harmonic transforms.

## Notes
* Intended to be a port of [ToyFold 2D](https://github.com/rhiju/toyfold2_rhiju/) to the 'real world', RNA folding in three dimensions
* Manipulation of SE(3) transforms for RNA inspired by [loop_close](https://github.com/rhiju/loop_close) repo (and `6D loop_close` implementation in Rosetta). 
* Starting with backbone conformational preferences drawn from ribosomal RNA

## Important TODO
* Dependence on frame choice -- see notes/ for possible solutions and cross-checks.
* KDE issue -- see notes/ for possible solutions and cross-checks.
