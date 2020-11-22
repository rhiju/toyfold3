function [C_eff,C_eff_err] = get_C_eff_overlap( step_types, TransformLibrary, NITER );
% [C_eff,C_eff_err] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
%
% Compute efficiently  C_eff based on overlap of
%  distributions going N/2 steps in forward and in reverse directions 
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
if ~exist('NITER','var') NITER = 1000; end;

i = floor(length(step_types)/2);

step_types_f = step_types(1:i);
all_pts_f = get_pts_forward( NITER, step_types_f, TransformLibrary);

step_types_r = step_types((i+1):end); % wait a minute...
pts = get_pts_forward( NITER, step_types_r, TransformLibrary);
all_pts_r = reverse_transform( pts );

[C_eff, C_eff_err] = get_C_eff_from_pts_6D(all_pts_f,all_pts_r);
