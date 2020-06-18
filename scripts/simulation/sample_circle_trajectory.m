function y = sample_circle_trajectory( step_types, TransformLibrary, all_pts_r );
% y = sample_circle_trajectory( step_types, TransformLibrary, all_pts_r );
%
% Sample forward trajectory, 'feeling' the effect of unbuilt links through
% the samples of reverse trajectories stored in all_pts_r
%
% INPUTS
%  step_types = list of steps ('BB',etc.) (Number of nucleotides N is length of this list plus 1 )
%  TransformLibrary = collection of TransformSets -- one must be BB.
%  all_pts_r = cell of Nmax arrays of reverse samples, each with 
%       [6 x NITER] points in 6D SE(3) space: x, y, z, v_x, v_y, v_z  
%       First, second, ... arrays correpond to transforms going from
%       *end* back 1 link, 2 links, etc.
%
% OUTPUT
%  y = struct with
%    t = [3 x N] coordinates of a random trace
%    R = [3 x 3 x N] orthonormal frames of a random trace
%


N = length(step_types);
x = zeros(3,N); m = zeros(3,3,N);
x(:,1) = [0,0,0]; % trajectory
m(:,:,1) = [1 0 0; 0 1 0; 0 0 1]; % orthonormal coordinate frame
for i = 1:N-1 % steps
    TransformSet = getfield(TransformLibrary,step_types{i});
    t = TransformSet.t;
    R = TransformSet.R;
    ntransforms = size(t,2);
    tic
    opts = struct();
    % positions going forward.
    for j = 1:ntransforms
        opts.t(:,j)  = x(:,i) + m(:,:,i)*t(:,j);
        opts.R(:,:,j)= m(:,:,i)*R(:,:,j);
    end
    opts = fill_T6_from_t_and_R( opts );
    pts_f = opts.T6;
    
    pts_r = all_pts_r{N-i}.T6; % reverse distribution
    
    s = get_kde_bandwidth( pts_r );
    p =  mvksdensity(pts_r,pts_f,'Bandwidth',s);
    
    p = p/sum(p);
    p_cumsum = cumsum( p+1e-10 );
    p_cumsum = p_cumsum / p_cumsum(end);
    
    q = rand(1);
    for idx = 1:ntransforms; if q <= p_cumsum(idx); break; end; end;
    x(:,i+1) = opts.t(:,idx);
    m(:,:,i+1) = opts.R(:,:,idx);
end
toc
x(:,N+1) = x(:,1);
m(:,:,N+1) = m(:,:,1);

clf
y = struct( 't',x,'R',m);
draw_trace(y,step_types)
hold on; 
plot3(0,0,0,'o');

