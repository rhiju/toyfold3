function draw_trace( t, R );
% draw_trace( y );
% draw_trace( t, R );
%
% Draw trace with coordinate frames.
%
% Input can be
%
%  t = [3 x N] coordinates of a  trace
%  R   = [3 x 3 x N] orthonormal frames of a  trace
%
% Or it can be
%
% y = struct with above fields.
%
% (C) R. Das, Stanford 2020

if isstruct( t );
    y = t;
    t = y.t;
    R = y.R;
end

plot3( t(1,:), t(2,:),t(3,:),'linew',2 ); hold on
for n = 1:size( t, 2 );
    xvec = [t(:,n), t(:,n) + 3*R(:,1,n)];
    yvec = [t(:,n), t(:,n) + 3*R(:,2,n)];
    zvec = [t(:,n), t(:,n) + 3*R(:,3,n)];
    plot3( zvec(1,:),zvec(2,:),zvec(3,:), 'k-','linew',2 );
    plot3( xvec(1,:),xvec(2,:),xvec(3,:), 'k-','linew',1 );
    plot3( yvec(1,:),yvec(2,:),yvec(3,:), 'k-','linew',0.5,'color',[0.5 0.5 0.5] );
end

axis equal
axis off
axis vis3d
set(gcf, 'PaperPositionMode','auto','color','white');
