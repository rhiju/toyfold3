function [C_eff_overlap_halfway,C_eff_overlap_halfway_error] = get_C_eff_overlap_halfway( N_overlap, all_pts_f, all_pts_r );
%
% Compute efficiently (but perhaps inaccurately?) C_eff based on overlap of
%  distributions going N/2 steps in forward and in reverse directions 
%
% Wrapper around get_C_eff_from_pts_6D.
%
%  N_overlap = lengths of circular RNA at which to evaluate overlap
%  all_pts_f = cell of Nmax arrays of forward samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%  all_pts_r = cell of Nmax arrays of reverse samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%
% (C) Rhiju Das, Stanford 2020
C_eff_overlap_halfway = [];
C_eff_overlap_halfway_error = [];
for k = 1:length(N_overlap)
    tic
    N = N_overlap(k);
    i = floor(N/2);
    [C_eff_overlap_halfway(k),C_eff_overlap_halfway_error(k)] = get_C_eff_from_pts_6D(all_pts_f{i},all_pts_r{N-i});
    toc
end