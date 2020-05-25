function [t,R] = get_transform( ctr1, M1, ctr2, M2);
% [t,R] = get_transform( ctr1, M1, ctr2, M2);
%
% Sanity check:
% ctr2 == ctr1 + M1*t
%
% Inputs
%  ctr1 = 3x1 x,y,z of input 
%  M1   = 3x3 input orthogonal frame
%  ctr2 = 3x1 x,y,z of output 
%  M2   = 3x3 output orthogonal frame
%
% Note format of M matrices:
%
%  M = 3x3 matrix such that M(:,1) = x unit vector, M(:,2) = y unit vector,
%        and M(:,3) = z unit vector
% Outputs
%
%  t = translation
%  R = rotation matrix
%
% (C) R. Das, Stanford 2020

t = inv(M1)*(ctr2 - ctr1);
R = inv(M1)*M2;

