function [t,R] = get_transform_library(ctr, M, chainbreak);
%
% get frame-to-frame transforms
%
% Inputs
%  ctr = [3 x N] coordinates of the trace
%  M   = [3 x 3 x N] orthonormal frames of the trace
%  chainbreak = [1 X N] is i to i+1 a chainbreak?
%
%  t = [3 x Nframes] library of translations from nt to nt. 
%  R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%
% (C) R. Das, Stanford 2020

t = []; % translations
R = []; % 3x3 rotation matrices
count = 0;
for n = find( ~chainbreak )
    % later need to put in a filter for chainbreaks!
    count = count+1;
    [t(:,count),R(:,:,count)] = get_transform( ctr(:,n), M(:,:,n), ctr(:,n+1), M(:,:,n+1));
end
