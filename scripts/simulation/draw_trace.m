function draw_trace( ctr, M );
% draw_trace( ctr, M );
%
% Draw trace with coordinate frames.
%
%  ctr = [3 x N] coordinates of a random trace
%  M   = [3 x 3 x N] orthonormal frames of a random trace
%
% (C) R. Das, Stanford 2020

plot3( ctr(1,:), ctr(2,:),ctr(3,:),'linew',2 ); hold on
for n = 1:size( ctr, 2 );
    xvec = [ctr(:,n), ctr(:,n) + 3*M(:,1,n)];
    yvec = [ctr(:,n), ctr(:,n) + 3*M(:,2,n)];
    zvec = [ctr(:,n), ctr(:,n) + 3*M(:,3,n)];
    plot3( zvec(1,:),zvec(2,:),zvec(3,:), 'k-','linew',2 );
    plot3( xvec(1,:),xvec(2,:),xvec(3,:), 'k-','linew',1 );
    plot3( yvec(1,:),yvec(2,:),yvec(3,:), 'k-','linew',0.5,'color',[0.5 0.5 0.5] );
end

axis equal
axis off
set(gcf, 'PaperPositionMode','auto','color','white');
