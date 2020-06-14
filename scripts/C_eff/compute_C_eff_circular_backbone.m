function [C_eff,C_eff_err] = compute_C_eff_circular_backbone(NITER, which_N, TransformLibrary);
% [C_eff,C_eff_err] = compute_C_eff_circular_backbone(NITER, which_N, TransformLibrary);
%
% Compute effective molarities for circularization for arbitrary
%  length N-mer, based on sampling few 1000 trajectories and
%  using multivariate 6D KDE estimate to get probability density
%  back at zero translation, no rotation.
% 
% Inputs
%  NITER = number of trajectories to sample for each length
%  which_N = what lengths of circles to compute C_eff for
%  TransformLibrary = collection of TransformSets -- one must be 'BB'.
%
% Outputs
%  C_eff = effective molarity (units of M) for each value of which_N
%  C_eff_err = effective molarity (units of M) for each value of which_N
%
% (C) R. Das, Stanford 2020

for i = 1:length(which_N)
    N = which_N(i); % number of steps
    tic
    step_types = repmat( {'BB'},N,1);
    [C_eff(i), C_eff_err(i)] = compute_C_eff_circular( NITER, step_types, TransformLibrary );
    toc
end


