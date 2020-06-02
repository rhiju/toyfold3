function C_eff = get_C_eff_from_pts_6D( pts_f, pts_r);
% C_eff = get_C_eff_from_pts_6D( pts_f, pts_r);
%
% INPUTS
%  pts_f = sampled (x,y,z,v_x,v_y_v_z) for forward chain segment
%  pts_r = sampled (x,y,z,v_x,v_y,v_z) for reverse chain segment
%
% OUTPUT
%  C_eff = Effective molarity for chain closure between forward samples
%            and reverse samples, calculated based on overlap of
%            distributions (over translations & rotations), and using
%            MATLAB's KDE estimator.
%
% (C) R. Das, Stanford University, 2020

isometric = 0; % leads to systematic underestimation
if isometric
    pts_f = convert_v_to_b(pts_f);
    pts_r = convert_v_to_b(pts_r);
end

s = get_kde_bandwidth( pts_f );

% phase space factor for integration with euler vectors
%  WAIT -- should I apply this or the inverse of it? 
%  Think through this more.

p = mvksdensity( pts_f, pts_r,'Bandwidth',s); 

if isometric
    C_eff = mean(p)/(1/(8*pi^2)*6.022e23/1e27 );
else
    % correct for phase space volume
    v = sqrt( sum(pts_r(:,4:6).^2, 2) ); % rotation angle
    w = (sin(v/2)./(v/2)).^2;
    C_eff = mean(p./(w/(8*pi^2)))/ (6.022e23/1e27);
end

function b = convert_v_to_b(v);
% isometric rescaling so uniform distribution in
% rotation vector space is uniform in sphere of radius
% (6*pi)^1/3
vnorm = sqrt( sum(v(:,4:6).^2,2) );
bnorm = ( 6 * (vnorm-sin(vnorm))).^(1/3);
b = v;
for n = 1:size( v,1); 
    b(n,4:6) = v(n,4:6) * b(n)/v(n); 
end;


