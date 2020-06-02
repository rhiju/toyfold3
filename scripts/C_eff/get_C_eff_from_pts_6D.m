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

s = get_kde_bandwidth( pts_f );

% phase space factor for integration with euler vectors
%  WAIT -- should I apply this or the inverse of it? 
%  Think through this more.

p = mvksdensity( pts_f, pts_r,'Bandwidth',s); 

v = sqrt( sum(pts_r(:,4:6).^2, 2) ); % rotation angle
w = (sin(v/2)./(v/2)).^2; 

%C_eff = mean(p)/(1/(8*pi^2)*6.022e23/1e27 );
C_eff = mean(p./(w/(8*pi^2)))/ (6.022e23/1e27);