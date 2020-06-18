function pts = fill_T6_from_t_and_R( pts );
% pts = fill_T6_from_t_and_R( pts );
%
% utility function to get from translations (3-vectors) and rotations (3x3)
%   to 6D vectors in SO(3) space.
%
%
% INPUT
% pts = struct with following fields
%    t = [3 x N] coordinates 
%    R = [3 x 3 x N] orthonormal frames 
% 
% OUTPUT
% pts = struct with t, R, and an additional field:
%    T6 = [6 x NITER] end-of-trajectory transforms in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%
%


% convert 3x3 rotation matrices to angle-axis (Euler vector)
pts_EV=SpinCalc('DCMtoEV',pts.R,0,0);  % output is unit vector, angle in degrees
pts_EV3 = [];
% convert to axis vector v_x,v_y,v_z; with length equal to rotation angle in radians
for n = 1:size( pts_EV,1); pts_EV3(n,:) = pts_EV(n,1:3) * pts_EV(n,4) * pi/180.0; end;

% points in 6D SE(3) space: x,y,z, v_x, v_y, v_z
pts.T6 = [pts.t', pts_EV3];
