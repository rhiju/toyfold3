function [C_eff_overlap_halfway,C_eff_overlap_halfway_error] = get_C_eff_overlap_halfway( N_overlap, all_pts_f, all_pts_r, N_overlap_offset );
%
% Compute efficiently  C_eff based on overlap of
%  distributions going N/2 steps in forward and in reverse directions 
%
% Be very careful in how you set up -- SCAN_LOOP_LENGTH provides a 
%  useful wrapper.
%
% Wrapper around get_C_eff_from_pts_6D.
%
% Inputs
%  N_overlap = lengths of circular RNA at which to evaluate overlap
%  all_pts_f = cell of Nmax arrays of forward samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%  all_pts_r = cell of Nmax arrays of reverse samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
% N_overlap_offset = index offset to apply when determining f/r to overlap
%                          (default 0)
%
% (C) Rhiju Das, Stanford 2020
if ~exist( 'N_overlap_offset' ) N_overlap_offset = 0; end;
C_eff_overlap_halfway = [];
C_eff_overlap_halfway_error = [];
for k = 1:length(N_overlap)
    tic
    N = N_overlap(k)          + N_overlap_offset;
    i = floor(N_overlap(k)/2) + N_overlap_offset;
    [C_eff_overlap_halfway(k),C_eff_overlap_halfway_error(k)] = get_C_eff_from_pts_6D(all_pts_f{i},all_pts_r{N-i});
    toc
end