function [x,m] = get_random_trace(N, t, R, drawit );

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

