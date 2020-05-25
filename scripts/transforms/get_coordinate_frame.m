function M = get_coordinate_frame( a, b, c );
%M = get_coordinate_frame( a, b, c );
%
% Get matrix with unit vectors for orthonormal frame defined
%  by a, b, c:
%
%  Point z from a to b
%  Point x perpendicular to z from a to c.
%  y is orthonormal to x,z.
%
% So:
%
%      a --  b
%             \
%              c
% leads to:
%
%   (y) is *into* page.
%     \  (z)
%      a --> b
%   (x)|      \
%      v       c
%
% INPUTS
%  a = 3-vector
%  b = 3-vector
%  c = 3-vector
%
% OUTPUT
%  M = 3x3 matrix such that M(:,1) = x unit vector, M(:,2) = y unit vector,
%        and M(:,3) = z unit vector
%
% (C) R. Das, Stanford 2020.

z = b - a;
z = z/norm(z);

x = c - b;
y = cross( z,x );
y = y/norm(y);

x = cross(y,z);

M = [x;y;z]';



