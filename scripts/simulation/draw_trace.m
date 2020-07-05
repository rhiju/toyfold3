function draw_trace( y, step_types );
% draw_trace( y, step_types );
%
% Draw trace with coordinate frames.
%
% Input should be
% y = struct with:
%  t = [3 x N] coordinates of a  trace
%  R = [3 x 3 x N] orthonormal frames of a  trace
%
% (C) R. Das, Stanford 2020

t = y.t;
R = y.R;

if ~exist( 'step_types', 'var' )
    % plot trace as usual.
    plot3( t(1,:), t(2,:),t(3,:),'linew',2 ); hold on
end

for n = 1:size( t, 2 );
    if exist( 'step_types', 'var' )
        % Draw trace, with colorcoding of types of transforms.
        if n < size( t, 2 )
            h = plot3( t(1,n+[0,1]), t(2,n+[0,1]),t(3,n+[0,1]),'linew',2 ); hold on
            colorcode = [0 0.4470 0.7410]; % blue/teal, default MATLAB color
            if strcmp(step_types{n},'Inline') colorcode = [1,0,0]; end;
            if strcmp(step_types{n},'BP') colorcode = [0.7,0.7,0.7]; end;
            set(h,'color',colorcode );
        end
    end
    
    % Draw coordinate frame
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
