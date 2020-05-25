function [t_rev, R_rev ] = reverse_transform( t, R);
% [t_rev, R_rev ] = reverse_transform( t, R);
% (C) R. Das, Stanford 2020
t_rev = inv(R)*(-t);
R_rev = inv(R);
