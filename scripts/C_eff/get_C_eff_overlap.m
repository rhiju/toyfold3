function [C_eff_f, C_eff_r,C_eff_f_err, C_eff_r_err] = get_C_eff_overlap( N, all_pts_f, all_pts_r, just_SO3 );
% [C_eff_f, C_eff_r,C_eff_f_err, C_eff_r_err] = get_C_eff_overlap( N, all_pts_f, all_pts_r );
%
%
% Inputs
%  N = length of circular RNA;
%  all_pts_f = cell of Nmax arrays of forward samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%  all_pts_r = cell of Nmax arrays of reverse samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%  just_SO3 = only look at rotation part of translation+rotation;
%
% Outputs
%
%  C_eff_f = effective molarity (units of M) for each value of #steps
%              forward from 1, 2, ... N
%  C_eff_r = effective molarity (units of M) for each value of #steps
%              reverse from 1, 2, ... N
%  C_eff_f_err = error in effective molarity (units of M) for each value of #steps
%              forward from 1, 2, ... N
%  C_eff_r_err = error in effective molarity (units of M) for each value of #steps
%              reverse from 1, 2, ... N
%
% (C) R. Das, Stanford University 2020
if ~exist( 'just_SO3','var') just_SO3 = 0; end;
for i = 1:N-1
    [C_eff_f(i),C_eff_f_err(i)] = get_C_eff_from_pts_6D(all_pts_f{i},all_pts_r{N-i},just_SO3);
    [C_eff_r(i),C_eff_r_err(i)] = get_C_eff_from_pts_6D(all_pts_r{i},all_pts_f{N-i},just_SO3);
end
[C_eff_f(N),C_eff_f_err(N)] = get_C_eff_from_pts_6D(all_pts_f{N},[0,0,0,0,0,0],just_SO3);
[C_eff_r(N),C_eff_r_err(N)] = get_C_eff_from_pts_6D(all_pts_r{N},[0,0,0,0,0,0],just_SO3);
