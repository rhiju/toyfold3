function [all_pts_f, all_pts_r] = get_all_pts( Nmax, NITER, t, R );
%
% Compute endpoint transforms for trajectories with step number 1 to Nmax.
% Note that number of 'nucleotides' at each step number N is N+1. We're
%  looking for that N+1-th nucleotide to return to first position.
%  
% Core data for 'everything' in C_eff and trajectory sampling calculations.
% 
% Inputs
%  Nmax  = max. length
%  NITER = number of trajectories to sample for each length
%  t = [3 x Nframes] library of translations from nt to nt. 
%  R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%
% Outputs
%  all_pts_f = cell of Nmax arrays of forward samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%  all_pts_r = cell of Nmax arrays of reverse samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%
% (C) R. Das, Stanford 2020
for i = 1:Nmax
    tic
    all_pts_f{i} = get_pts_forward( NITER, i+1, t, R);
    % don't really need to repeat above -- could just reverse transforms
    all_pts_r{i} = get_pts_reverse( NITER, i+1, t, R); 
    toc
end

