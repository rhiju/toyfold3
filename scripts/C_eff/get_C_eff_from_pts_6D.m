function [C_eff,C_eff_error,C_eff_samples] = get_C_eff_from_pts_6D( pts_f, pts_r, just_SO3, both_dirs);
% [C_eff,C_eff_error] = get_C_eff_from_pts_6D( pts_f, pts_r, just_SO3);
%
% INPUTS
%  pts_f = sampled (x,y,z,v_x,v_y_v_z) for forward chain segment (could be
%              6D tensor, OR struct holding tensor in T6 field.)
%  pts_r = sampled (x,y,z,v_x,v_y,v_z) for reverse chain segment (could be
%              6D tensor, OR struct holding tensor in T6 field.)
%  just_SO3 = only look at rotation part of translation+rotation (default
%  0)
%  both_dirs = do KDE with forward/reverse *and* reverse/forward to
%                get better estimate of C_eff and its error. (default 0,
%                but highly recommended)
%
% OUTPUT
%  C_eff = Effective molarity for chain closure between forward samples
%            and reverse samples, calculated based on overlap of
%            distributions (over translations & rotations), and using
%            MATLAB's KDE estimator.
%  C_eff_err = error estimate in C_eff. 
%                NOTE: by default, only includes statistical error from
%                   sampling points in pts_r. 
%                NOTE: If both_dirs=1, stat. error from sampling pts_f is  
%                   also included, as well as systematic error estimated
%                   from deviation between forward/reverse and reverse/
%                   forward calculations.
%  C_eff_samples = samples used to estimate C_eff [same number of values as pts_r]
%                    NOTE: C_eff = mean( C_eff_samples). 
%                    NOTE: returns empty if both_dirs = 1.
%
% (C) R. Das, Stanford University, 2020
if ~exist( 'just_SO3','var') just_SO3 = 0; end;
if ~exist( 'both_dirs','var') both_dirs = 0; end;

if both_dirs
    [C_eff_f, C_eff_err_f] = get_C_eff_from_pts_6D(pts_f,pts_r);
    [C_eff_r, C_eff_err_r] = get_C_eff_from_pts_6D(pts_r,pts_f);
    C_eff = sqrt(C_eff_f * C_eff_r);
    
    % statistical error
    C_eff_relerr_f = C_eff_err_f/C_eff_f;
    C_eff_relerr_r = C_eff_err_r/C_eff_r;
    % systematic error
    C_eff_diff = 0.5 * abs(log(C_eff_f/C_eff_r));
    
    C_eff_relerr = sqrt(C_eff_relerr_f^2 + C_eff_relerr_r^2 + C_eff_diff^2);
    C_eff_error = C_eff_relerr * C_eff;
    C_eff_samples = [];
    return
end

if isstruct( pts_f ) 
    if ~isfield(pts_f,'T6'); pts_f = fill_T6_from_t_and_R(pts_f); end;
    pts_f = pts_f.T6; 
end;
if isstruct( pts_r ) 
    if ~isfield(pts_r,'T6'); pts_r = fill_T6_from_t_and_R(pts_r); end;
    pts_r = pts_r.T6; 
end;
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
    C_eff_samples = p/(1/(8*pi^2)*6.022e23/1e27 );
else
    % correct for phase space volume
    v = sqrt( sum(pts_r(:,end-2:end).^2, 2) ); % rotation angle
    w = (sin(v/2)./(v/2)).^2;
    w( find( v == 0 ) ) = 1;
    C_eff_samples = p./(w/(8*pi^2))/ (6.022e23/1e27);
end

if just_SO3
    % don't need to conversion to molarity.
    C_eff_samples = C_eff_samples * (6.022e23/1e27);
end
    
C_eff = mean(C_eff_samples);
C_eff_error = std(C_eff_samples)/sqrt(length(p));


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


