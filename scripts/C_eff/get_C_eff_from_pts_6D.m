function [C_eff,C_eff_error] = get_C_eff_from_pts_6D( pts_f, pts_r, just_SO3);
% [C_eff,C_eff_error] = get_C_eff_from_pts_6D( pts_f, pts_r, just_SO3);
%
% INPUTS
%  pts_f = sampled (x,y,z,v_x,v_y_v_z) for forward chain segment (could be
%              6D tensor, OR struct holding tensor in T6 field.)
%  pts_r = sampled (x,y,z,v_x,v_y,v_z) for reverse chain segment (could be
%              6D tensor, OR struct holding tensor in T6 field.)
%  just_SO3 = only look at rotation part of translation+rotation;
%
% OUTPUT
%  C_eff = Effective molarity for chain closure between forward samples
%            and reverse samples, calculated based on overlap of
%            distributions (over translations & rotations), and using
%            MATLAB's KDE estimator.
%  C_eff_err = error estimate in C_eff. NOTE: only includes error from
%                   sampling points in pts_r -- probably should  bootstrap 
%                   over pts_f too?
%
% (C) R. Das, Stanford University, 2020
if ~exist( 'just_SO3','var') just_SO3 = 0; end;
if isstruct( pts_f ) pts_f = pts_f.T6; end;
if isstruct( pts_r ) pts_r = pts_r.T6; end;
C_eff_error = NaN;
    
isometric = 0; % leads to systematic underestimation
if isometric
    pts_f = convert_v_to_b(pts_f);
    pts_r = convert_v_to_b(pts_r);
end
if just_SO3
    pts_f = pts_f(:,4:6);
    pts_r = pts_r(:,4:6);
end

s = get_kde_bandwidth( pts_f );
p = mvksdensity( pts_f, pts_r,'Bandwidth',s);

if isometric
    C_eff = mean(p)/(1/(8*pi^2)*6.022e23/1e27 );
    C_eff_error = std(p)/(1/(8*pi^2)*6.022e23/1e27 )/sqrt(length(p));
else
    % correct for phase space volume
    v = sqrt( sum(pts_r(:,end-2:end).^2, 2) ); % rotation angle
    w = (sin(v/2)./(v/2)).^2;
    w( find( v == 0 ) ) = 1;
    C_eff = mean(p./(w/(8*pi^2)))/ (6.022e23/1e27);
    C_eff_error = std(p./w)/(1/(8*pi^2)*6.022e23/1e27 )/sqrt(length(p));
end

if just_SO3
    % don't need to conversion to molarity.
    C_eff = C_eff * (6.022e23/1e27);
    C_eff_error = C_eff_error * (6.022e23/1e27);
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


