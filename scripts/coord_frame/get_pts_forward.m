function [pts_x,pts_m] = get_pts_forward( NITER, N, t, R);

for i = 1:NITER
    [x,m] = get_random_trace(N, t, R, 0);
    pts_x(:,i) = x(:,N);
    pts_m(:,:,i) = m(:,:,N);
end