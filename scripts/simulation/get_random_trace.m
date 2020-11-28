function y = get_random_trace(step_types, TransformLibrary, drawit );
% y = get_random_trace(step_types, TransformLibrary,  drawit );
%
% Coordinates of random trace,
%  starting from (0,0,0) and standard rotation frame.
%
% Inputs
%  step_types = list of steps ('BB',etc.) (Number of nucleotides N is length of this list plus 1 )
%  TransformLibrary = collection of TransformSets with field names corresponding to
%       step_types.
%  drawit = draw random trace with coordinate frames [default = 1]
%
% Output
%  y = struct with
%    t = [3 x N] coordinates of a random trace
%    R = [3 x 3 x N] orthonormal frames of a random trace
%
% (C) R. Das, Stanford 2020

if ~exist( 'drawit', 'var' ) drawit = 1; end;
N = length( step_types ) + 1;
x = zeros(3,N); m = zeros(3,3,N);
x(:,1) = [0,0,0]; % trajectory
m(:,:,1) = [1 0 0; 0 1 0; 0 0 1]; % orthonormal coordinate frame
for n = 2:N
    transforms = getfield( TransformLibrary, step_types{n-1} );
    t = transforms.t;
    R = transforms.R;
    ntransforms = size(t,2);
    j = randi(ntransforms);
    [x(:,n), m(:,:,n)] = apply_transform( x(:,n-1), m(:,:,n-1), t(:,j), R(:,:,j) );
    x(:,n)  = x(:,n-1) + m(:,:,n-1)*t(:,j);
    m(:,:,n)= m(:,:,n-1)*R(:,:,j);
end
y = struct( 't',x,'R',m);
if drawit
    cla
    draw_trace(y,step_types)
end




