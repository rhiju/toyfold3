function s = get_kde_bandwidth( pts );
% get_kde_bandwidth( pts );
% 
% Using Silverman't rule of thumb, as 
%  described in MATLAB's mvskdensity reference page!
%
% (C) R. Das, Stanford University 2019
%
% [Same function as in Toyfold2 repo]

d = size(pts,1);
n = size(pts,2);
s = std( pts )*(4/(n+2)/d)^(1/(n+4)); % bandwidth estimator