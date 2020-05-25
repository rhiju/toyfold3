function [t,R] = get_transform_library(ctr, M);
%
% get frame-to-frame transforms
%
% Inputs
%  ctr = [3 x N] coordinates of the trace
%  M   = [3 x 3 x N] orthonormal frames of the trace
%
%  t = [3 x Nframes] library of translations from nt to nt. 
%  R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%
% (C) R. Das, Stanford 2020

t = []; % translations
R = []; % 3x3 rotation matrices
for n = 1:(size( ctr, 2)-1);
    % later need to put in a filter for chainbreaks!
    [t(:,n),R(:,:,n)] = get_transform( ctr(:,n), M(:,:,n), ctr(:,n+1), M(:,:,n+1));
end
