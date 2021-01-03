function y = sample_circle_trajectories( step_types, TransformLibrary, NTRAJ, NITER );
% y = sample_circle_trajectories( step_types, TransformLibrary, NTRAJ, NITER );
%
% Sample forward trajectory, 'feeling' the effect of unbuilt links through
% the samples of reverse trajectories stored in all_pts_r
%
% INPUTS
%  step_types = list of steps ('BB',etc.) (Number of nucleotides N is length of this list plus 1 )
%  TransformLibrary = collection of TransformSets -- one must be BB.
%  NTRAJ = number of trajectories (default 1)
%  NITER = number of iterations for computing 6D convolutions(default 1000)
%
% OUTPUT
%  y = cell of structs with
%    t = [3 x N] coordinates of a random trace
%    R = [3 x 3 x N] orthonormal frames of a random trace
%
% (C) R. Das, Stanford University 2021

if ~exist( 'NTRAJ','var') NTRAJ = 1; end;
if ~exist( 'NITER','var') NITER = 1000; end;

all_pts_r = get_all_pts_r(step_types, TransformLibrary, NITER);

y = {};
for n = 1:NTRAJ
    fprintf( 'Doing trajectory %d of %d...\n',n.NTRAJ);
    y{n} = sample_circle_trajectory( step_types, TransformLibrary, all_pts_r, 0 );
end


function  all_pts_r = get_all_pts_r(step_types, TransformLibrary, NITER);
% stolen from get_all_pts()
for i = 1:length(step_types)
    tic
    step_types_r = step_types(end-i+1:end);
    pts = get_pts_forward( NITER, step_types_r, TransformLibrary);
    all_pts_r{i} = reverse_transform( pts ); 
    toc
end