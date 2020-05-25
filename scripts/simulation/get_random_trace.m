function [x,m] = get_random_trace(N, t, R, drawit );
% [x,m] = get_random_trace(N, t, R, drawit );
%
% Coordinates of random trace,
%  starting from (0,0,0) and standard rotation frame.
%
% Inputs
%  N = number of nucleotides. (Number of steps is this number -1 )
%  t = [3 x Nframes] library of translations from nt to nt. 
%  R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%  drawit = draw random trace with coordinate frames [default = 1]
%
% Output
%  x = [3 x N] coordinates of a random trace
%  m = [3 x 3 x N] orthonormal frames of a random trace
%
% (C) R. Das, Stanford 2020

if ~exist( 'drawit', 'var' ) drawit = 1; end;

x = zeros(3,N); m = zeros(3,3,N);
x(:,1) = [0,0,0]; % trajectory
m(:,:,1) = [1 0 0; 0 1 0; 0 0 1]; % orthonormal coordinate frame
ntransforms = size(t,2);
for n = 2:N
    j = randi(ntransforms);
    x(:,n)  = x(:,n-1) + m(:,:,n-1)*t(:,j);
    m(:,:,n)= m(:,:,n-1)*R(:,:,j);
end

if drawit
    cla
    draw_trace(x,m)
end

