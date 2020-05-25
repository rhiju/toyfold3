function pts = get_pts_forward( NITER, N, t, R);
%
% NOTE -- this is redundant with get_pts_reverse --
%  just need to write a function to reverse transforms that takes
%  the 6D vector below (rather than trans, 3x3 rotation).
%
% Randomly sample a bunch of trajectories and record the reverse transform
% to get from the endpoints to the start
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
% Inputs
%  NITER = draw random trace with coordinate frames [default = 1]
%  N = number of nucleotides. (Number of steps is this number -1 )
%  t = [3 x Nframes] library of translations from nt to nt. 
%  R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%
% Outputs
%  pts = [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%
% (C) R. Das, Stanford 2020
for i = 1:NITER
    [x,m] = get_random_trace(N, t, R, 0);
    [x_rev,m_rev] = reverse_transform( x(:,N), m(:,:,N) );
    pts_x(:,i) = x_rev;
    pts_m(:,:,i) = m_rev;
end

% convert 3x3 rotation matrices to angle-axis (Euler vector)
pts_EV=SpinCalc('DCMtoEV',pts_m,0,0);  % output is unit vector, angle in degrees
pts_EV3 = [];
% convert to axis vector v_x,v_y,v_z; with length equal to rotation angle in radians
for n = 1:size( pts_EV,1); pts_EV3(n,:) = pts_EV(n,1:3) * pts_EV(n,4) * pi/180.0; end;

% points in 6D SE(3) space: x,y,z, v_x, v_y, v_z
pts = [pts_x', pts_EV3];
