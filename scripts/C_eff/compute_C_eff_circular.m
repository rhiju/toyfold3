function C_eff = compute_C_eff_circular(NITER, which_N, t, R);
% C_eff = compute_C_eff_circular(NITER, which_N, t, R);
%
% Compute effective molarities for circularization for arbitrary
%  length N-mer, based on sampling few 1000 trajectories and
%  using multivariate 6D KDE estimate to get probability density
%  back at zero translation, no rotation.
% 
% Inputs
%  NITER = number of trajectories to sample for each length
%  which_N = what lengths of circles to compute C_eff for
%  t = [3 x Nframes] library of translations from nt to nt. 
%  R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%
% Outputs
%  C_eff = effective molarity (units of M) for each value of which_N
%
% (C) R. Das, Stanford 2020

for i = 1:length(which_N)
    N = which_N(i); % number of steps
    tic
    pts_f = get_pts_forward( NITER, N+1, t, R);
    toc
    %cla;plot3( pts_f(1,:),pts_f(2,:),pts_f(3,:),'o'); axis equal
    % then collect histograms
        
    s = get_kde_bandwidth( pts_f );
    pts_r = [0,0,0,0,0,0];

    p(i) =  mvksdensity(pts_f,pts_r,'Bandwidth',s)';
end

C_eff = p/(1/(8*pi^2)*6.022e23/1e27 );
