%load toyfold3_test.mat

N = 12; % steps.

x = zeros(3,N); m = zeros(3,3,N);
x(:,1) = [0,0,0]; % trajectory
m(:,:,1) = [1 0 0; 0 1 0; 0 0 1]; % orthonormal coordinate frame
ntransforms = size(t,2);
tic
for i = 1:N-1 % steps
    % positions going forward.
    for j = 1:ntransforms
        opts_x(:,j)  = x(:,i) + m(:,:,i)*t(:,j);
        opts_m(:,:,j)= m(:,:,i)*R(:,:,j);
    end
    % convert 3x3 rotation matrices to angle-axis (Euler vector)
    opts_EV=SpinCalc('DCMtoEV',opts_m,0,0);  % output is unit vector, angle in degrees
    opts_EV3 = [];
    % convert to axis vector v_x,v_y,v_z; with length equal to rotation angle in radians
    for n = 1:size( opts_EV,1); opts_EV3(n,:) = opts_EV(n,1:3) * opts_EV(n,4) * pi/180.0; end;
    pts_f = [opts_x', opts_EV3];
    
    pts_r = all_pts_r{N-i}; % reverse distribution
    
    s = get_kde_bandwidth( pts_r );
    p =  mvksdensity(pts_r,pts_f,'Bandwidth',s);
    
    p = p/sum(p);
    p_cumsum = cumsum( p+1e-10 );
    p_cumsum = p_cumsum / p_cumsum(end);
    
    q = rand(1);
    for idx = 1:ntransforms; if q <= p_cumsum(idx); break; end; end;
    x(:,i+1) = opts_x(:,idx);
    m(:,:,i+1) = opts_m(:,:,idx);
end
toc
x(:,N+1) = x(:,1);
m(:,:,N+1) = m(:,:,1);

clf
draw_trace(x,m)
hold on; 
plot3(0,0,0,'o');

