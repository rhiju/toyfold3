function out = scan_loop_length( loop_lengths, NITER, TransformLibrary, motif );
% out = scan_loop_length( loop_lengths, NITER, motif );
%
% Scan length of a contiguous loop incorporated into a motif.
% 
% Compute efficiently C_eff based on overlap of
%    distributions going N/2 steps in forward and in reverse directions.
%
% Wrapper around GET_C_EFF_OVERLAP_HALFWAY
%
% Inputs
%  loop_lengths = lengths of loops to try, e.g. 1:20. Values must be greater than 1.
%  NITER = number of trajectories to sample for each length
%  TransformLibrary = collection of TransformSets -- one must be BB.
%  motif = [optional] step_types, e.g. {'BP','BB'} to seed a hairpin.
%
% Output
%  out = struct with fields:
%   out.C_eff: circularization C_eff with motif extended by each N in loop_lengths
%   out.C_eff_err:  approximate error in above.
%
% (C) R. Das, Stanford University, 2020

if ~exist( 'motif','var') motif = {}; end;

N_overlap_offset = length( motif );
step_types = [motif,repmat({'BB'},1,max(loop_lengths))];

[all_pts_f, all_pts_r] = get_all_pts( step_types, NITER, TransformLibrary );

[C_eff_overlap_halfway,C_eff_overlap_halfway_error] = get_C_eff_overlap_halfway_wrapper( loop_lengths, all_pts_f, all_pts_r, N_overlap_offset );

out.C_eff = C_eff_overlap_halfway;
out.C_eff_err = C_eff_overlap_halfway_error;



