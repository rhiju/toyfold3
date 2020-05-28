# Toyfold 3
Calculations and simulations of RNA tertiary structure in 3D

_Author: Rhiju Das, 2020_

## Goals
* Allow rapid simulations to test partitition function calcs -- including tertiary structure -- in MATLAB.
* Actual predictions for `C_eff` and structural ensembles for circular RNA
* Actual predictions for nearest neighbor energies and structural ensembles for helices
* Test factorization of free energy in terms of local motif energy (including K_d); motif modularity costs; tertiary closure.
* Test simple analytic representations of translation/rotation SE(3) distributions in tertiary closure costs, including Gaussian models and harmonic transforms.

![Example trace of random conformation](notes/ToyFold3D_NOTES.rtfd/Screen%20Shot%202020-05-25%20at%2012.13.29%20PM.png)

## Notes
* Intended to be a port of [ToyFold 2D](https://github.com/rhiju/toyfold2_rhiju/) to the 'real world', RNA folding in three dimensions
* Manipulation of SE(3) transforms for RNA inspired by [loop_close](https://github.com/rhiju/loop_close) repo (and `6D loop_close` implementation in Rosetta). 
* Starting with backbone conformational preferences drawn from ribosomal 23S RNA (4YBB).

## Getting Started
* Add to your paths: `scripts/` and subdirectories.
* Check `notes/` for TextEdit viewable notes on current modeling status
* Check `runs/` for example commands. First load data saved via `load toyfold3_test.mat`. Then check the commands in, say, `sample_circle_trajectory_test.m`. Setup of the needed variables are in `toyfold3_test.m`. 

![Example of a circular RNA trace](notes/ToyFold3D_NOTES.rtfd/Screen%20Shot%202020-05-25%20at%203.26.51%20PM.png)


## TODO
* KDE issue -- see notes/ for possible solutions and cross-checks. Probably need to work out *grid* base sampling; can see Toyfold 2d [get\_sample\_trajectories\_GRID.m](https://github.com/rhiju/toyfold2_rhiju/blob/master/scripts/sample_trajectories/get_sample_trajectories_GRID.m) for example.
* Dependence on frame choice -- see `notes/` for possible solutions and cross-checks.
