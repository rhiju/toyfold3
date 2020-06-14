function C_eff = compute_C_eff_circular(NITER, which_N, TransformLibrary);
% C_eff = compute_C_eff_circular(NITER, which_N, TransformLibrary);
%
% Compute effective molarities for circularization for arbitrary
%  length N-mer, based on sampling few 1000 trajectories and
%  using multivariate 6D KDE estimate to get probability density
%  back at zero translation, no rotation.
% 
% Inputs
%  NITER = number of trajectories to sample for each length
%  which_N = what lengths of circles to compute C_eff for
%  TransformLibrary = collection of TransformSets -- one must be BB.
%
% Outputs
%  C_eff = effective molarity (units of M) for each value of which_N
%
% (C) R. Das, Stanford 2020

for i = 1:length(which_N)
    N = which_N(i); % number of steps
    tic
    step_types = repmat( {'BB'},N,1);
    pts_f = get_pts_forward( NITER, step_types, TransformLibrary);
    toc
    %cla;plot3( pts_f.t(1,:),pts_f.t(2,:),pts_f.t(3,:),'o'); axis equal
    % then collect histograms
    pts_f = pts_f.T6; % 6D tensor
    
    s = get_kde_bandwidth( pts_f );
    pts_r = [0,0,0,0,0,0];

    % note there is no sinc(v/2)^2 volume element factor at (0,0,0) in
    % SO(3)!
    tic
    p(i) =  mvksdensity(pts_f,pts_r,'Bandwidth',s)';
    toc
end

C_eff = p/(1/(8*pi^2)*6.022e23/1e27 );
