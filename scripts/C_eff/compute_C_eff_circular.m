function [C_eff, C_eff_err] = compute_C_eff_circular( NITER, step_types, TransformLibrary, NBOOT );
% [C_eff,C_eff_err] = compute_C_eff_circular(NITER, which_N, TransformLibrary);
%
% Compute effective molarities for circularization for trajectory, based 
%  on sampling few 1000 trajectories and using multivariate 6D KDE estimate 
%  to get probability density
%  back at zero translation, no rotation.
% 
% Inputs
%  NITER = number of trajectories to sample for each length
%  step_types = list of steps ('BB',etc.) (Number of nucleotides N is length of this list plus 1 )
%  TransformLibrary = collection of TransformSets -- one must be 'BB'.
%
% Outputs
%  C_eff = effective molarity (units of M) for each value of which_N
%  C_eff_err = effective molarity (units of M) for each value of which_N
%
% (C) R. Das, Stanford 2020

pts_f = get_pts_forward( NITER, step_types, TransformLibrary);

%cla;plot3( pts_f.t(1,:),pts_f.t(2,:),pts_f.t(3,:),'o'); axis equal
%hold on; plot3(0,0,0,'ko','markersize',10); pause;
% then collect histograms
pts_f = pts_f.T6; % 6D tensor
pts_r = [0,0,0,0,0,0];

tic
Nclusters = 1;
if Nclusters > 0
    C_eff = get_C_eff_clustered( pts_f, pts_r, Nclusters );
else
    s = get_kde_bandwidth( pts_f);
    % note there is no sinc(v/2)^2 volume element factor at (0,0,0) in SO(3)!
    p = mvksdensity(pts_f,pts_r,'Bandwidth',s)';
    C_eff = p/(1/(8*pi^2)*6.022e23/1e27 );
end
toc
 
% 6D Gaussian approximation near origin...
% gp = find( abs(pts_f(:,1))<0.5 & abs(pts_f(:,2))<0.5 & abs(pts_f(:,3))<0.5 );
% c = cov( pts_f(gp,:) );
% C_eff_gaussian = 1/(2*pi)^3/sqrt(det(c))/(1/(8*pi^2)*6.022e23/1e27)


C_eff_err = NaN;
% bootstrap?
%Nsamples = length( pts_f );
%randi( size( pts_f,1) )



%%%%%%%%%%%%%
% experimental...
function  C_eff = get_C_eff_clustered( pts_f, pts_r, Nclusters );

cidx = kmeans( pts_f, Nclusters );
C_eff = 0;
for q = 1:Nclusters
    s = get_kde_bandwidth( pts_f(cidx==q,:) );
    % note there is no sinc(v/2)^2 volume element factor at (0,0,0) in
    % SO(3)!
    p = mvksdensity(pts_f(cidx==q,:),pts_r,'Bandwidth',s)';
    C_eff = C_eff + p/(1/(8*pi^2)*6.022e23/1e27 );
end
