function [C_eff,C_eff_err] = get_C_eff_overlap( step_types, TransformLibrary, NITER );
% [C_eff,C_eff_err] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
%
% Compute efficiently  C_eff based on overlap of
%  distributions going N/2 steps in forward and in reverse directions 
%
% Inputs
%  step_types = {'BP','BB','BB','BB','BB','BB'} is tetraloop.
%  TransformLibrary = collection of TransformSets -- one must be BB.
%  NITER = number of trajectories to sample for each length
%
% Outputs 
%  C_eff = effective molarity (units of M) for overlap
%  C_eff_err = error in effective molarity (units of M) for overlap
%
% (C) R. Das, Stanford, 2020

if ~exist('NITER','var') NITER = 1000; end;

i = floor(length(step_types)/2);

step_types_f = step_types(1:i);
all_pts_f = get_pts_forward( NITER, step_types_f, TransformLibrary);

step_types_r = step_types((i+1):end); % wait a minute...
pts = get_pts_forward( NITER, step_types_r, TransformLibrary);
all_pts_r = reverse_transform( pts );

%[C_eff, C_eff_err] = get_C_eff_from_pts_6D(all_pts_f,all_pts_r);

[C_eff_f, C_eff_err_f] = get_C_eff_from_pts_6D(all_pts_f,all_pts_r);
[C_eff_r, C_eff_err_r] = get_C_eff_from_pts_6D(all_pts_r,all_pts_f);
C_eff = sqrt(C_eff_f * C_eff_r);

% statistical error
C_eff_relerr_f = C_eff_err_f/C_eff_f;
C_eff_relerr_r = C_eff_err_r/C_eff_r;
% systematic error
C_eff_diff = 0.5 * abs(log(C_eff_f/C_eff_r));

C_eff_relerr = sqrt(C_eff_relerr_f^2 + C_eff_relerr_r^2 + C_eff_diff^2);
C_eff_err = C_eff_relerr * C_eff;


