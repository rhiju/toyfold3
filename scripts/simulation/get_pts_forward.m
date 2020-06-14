function pts = get_pts_forward( NITER, step_types, TransformLibrary);
% pts = get_pts_forward( NITER, step_types, TransformLibrary);
%
% Randomly sample a bunch of trajectories and record their
%  endpoints. Using 
%
%     x, y, z, v_x, v_y, v_z
%
%  where x, y, and z give translation in Angstroms, and 
%  v_x, v_y, and v_z give axis vector of the rotation whose length 
%    gives rotation (in radians). 
%
% To convert the rotation to other formats like
%    Euler angles or 3x3 rotation matrix, use SpinCalc function!
%
%
% Inputs
%  NITER = number of trajectories to sample for each length
%  step_types = list of steps ('BB',etc.) (Number of nucleotides N is length of this list plus 1 )
%  TransformLibrary = collection of TransformSets -- one must be BB.

% Outputs
%  pts = struct with 3 useful fields:
%    t = [3 x NITER] coordinates of a random trace
%    R = [3 x 3 x NITER] orthonormal frames of a random trace
%    T6 = [6 x NITER] end-of-trajectory transforms in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%
% (C) R. Das, Stanford 2020

N = length( step_types )  + 1;
for i = 1:NITER
    y = get_random_trace(step_types, TransformLibrary, 0);
    pts.t(:,i) = y.t(:,N);
    pts.R(:,:,i) = y.R(:,:,N);
end

pts = fill_T6_from_t_and_R( pts );

