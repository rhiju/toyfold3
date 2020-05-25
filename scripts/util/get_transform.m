function [t,R] = get_transform( ctr1, M1, ctr2, M2);
% Outputs
%  t = translation
%  R = rotation matrix

t = inv(M1)*(ctr2 - ctr1);
R = inv(M1)*M2;

% Sanity check:
% ctr2 == ctr1 + M1*t
