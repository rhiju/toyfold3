function [t_rev, R_rev] = reverse_transform(t, R);
% y_rev = reverse_transform( y );
% [t_rev, R_rev] = reverse_transform(t, R);
%
%  t = [3 x N] coordinates 
%  R = [3 x 3 x N] orthonormal frames 
% 
%
% (C) R. Das, Stanford 2020
if isstruct( t )
    y = t;
    [t_rev, R_rev] = reverse_transform( y.t, y.R);
    y_rev = struct();
    y_rev.t = t_rev;
    y_rev.R = R_rev;
    if isfield( y, 'T6' );
        y_rev = fill_T6_from_t_and_R( y_rev );
    end
    t_rev = y_rev;
    clear R_rev;
    return;
end

if length( size( R ) )>2
    for i = 1:size( R, 3 )
        [t_rev(:,i), R_rev(:,:,i)] = reverse_transform( t(:,i), R(:,:,i) );
    end
    return;
end

t_rev = inv(R)*(-t);
R_rev = inv(R);
