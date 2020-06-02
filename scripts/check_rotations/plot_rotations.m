function v = plot_rotations( phi, theta, psi, rotation_vector, gp, use_axis_angle );
% v = plot_rotations( phi, theta, psi, rotation_vector_radians, gp, use_axis_angle );
%
%  INPUT:
%   phi   = euler angle in degrees  (first rotation, Z)
%   theta = euler angle in degrees  (second rotation, X)
%   psi   = euler angle in degrees  (third rotation, Z)
%   rotation_vector = rotation_vector, magnitude is angle of rotation in
%                        [0,180]  (i.e., in degrees)
%   gp    = 'good points', indices of points to plot
%   use_axis_angle = plot with axis_angle (in degrees), not
%                       euler
%
%  OUTPUT:
%   v = rotation_vector, if computed. (in *degrees*)
%

if nargin <= 3
    D = phi;
    phi = D(:,1);
    theta = D(:,2);
    psi = D(:,3);
    if nargin > 1; 
        gp = theta;
    end
end
if exist( 'gp', 'var' ); gp == 0 | gp == 1; axis_angle = gp; clear gp; end;
if ~exist( 'gp', 'var' ); gp = [1:length(phi)]; end;
if ~exist( 'use_axis_angle', 'var' ) use_axis_angle = 0; end;
v = [];    
if ~use_axis_angle
   plot3( phi(gp),cos(theta(gp)*pi/180),psi(gp),'.');
   xlabel( '\alpha' ); ylabel( 'cos \beta' ); zlabel( '\gamma' );
else
   if ~isempty( rotation_vector )
    v = rotation_vector;
   else
    v = get_rotation_vector_from_euler( phi,theta, psi );
   end
   plot3( v(gp,1), v(gp,2), v(gp,3),'.');
   xlabel( 'v_x' ); ylabel( 'v_y' ); zlabel( 'v_z' );
end
set(gcf, 'PaperPositionMode','auto','color','white');
axis vis3d

