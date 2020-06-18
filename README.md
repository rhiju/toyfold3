# Toyfold 3
Calculations of RNA tertiary structure ensembles and free energies in 3D

_Author: Rhiju Das, 2020_

## Goals
* Allow rapid simulations to test partitition function calcs -- including tertiary structure -- in MATLAB.
* Actual predictions for `C_eff` and structural ensembles for circular RNA
* Explore sampling of ribose rings as a 'sub problem' in atomic reconstruction.
* Actual predictions for nearest neighbor energies and structural ensembles for helices
* Test factorization of free energy in terms of local motif energy (including K_d); motif modularity costs; tertiary closure.
* Test simple analytic representations of translation/rotation SE(3) distributions in tertiary closure costs, including Gaussian models and harmonic transforms.

![Example trace of random conformation](notes/01_ToyFold3D_NOTES.rtfd/Screen%20Shot%202020-05-25%20at%2012.13.29%20PM.png)

## Notes
* Intended to be a port of [ToyFold 2D](https://github.com/rhiju/toyfold2_rhiju/) to the 'real world', RNA folding in three dimensions
* Manipulation of SE(3) transforms for RNA inspired by [loop_close](https://github.com/rhiju/loop_close) repo (and `6D loop_close` implementation in Rosetta). 
* Starting with backbone conformational preferences drawn from ribosomal 23S RNA (4YBB).

## Getting Started
* Add to your paths: `scripts/` and subdirectories.
* Check `notes/` for TextEdit viewable notes on current modeling status
* Check `runs/` for example commands. First load data saved via `load toyfold3_test.mat`. Then check the commands in, say, `toyfold3_test.m`, particularly the code block at the end that runs `sample_circle_trajectory`. 

![Example of a circular RNA trace](notes/01_ToyFold3D_NOTES.rtfd/Screen%20Shot%202020-05-25%20at%203.26.51%20PM.png)


