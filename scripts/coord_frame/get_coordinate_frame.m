function M = get_coordinate_frame( a, b, c );
% Point z from a to b
% Point x perpendicular to z from a to c.
% y is orthonormal.
%
%         z
%      a --> b
%     x|      \
%      v       c
z = b - a;
z = z/norm(z);

x = c - b;
y = cross( z,x );
y = y/norm(y);

x = cross(y,z);

M = [x;y;z]';



