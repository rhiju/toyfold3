d = pdbread( 'data/4ybb_DIII.pdb');

%%
% grab coordinate frames O5'
%
%    O5' <-- center
%     z\      z points from O5' to C5'
%       C5'   x points perpendicular to z towards C4'
%       |    [y completes orthonormal frame.]
%       C4'
%
count = 0;
resnum = [];
chain = '';
o5prime = []; c5prime = []; c4prime = [];
xyz = []; % collect all xyz, just for checks
for i = 1:length( d.Model.Atom )
    atom = d.Model.Atom(i);
    if strcmp( atom.AtomName, 'C5''' )
        count = count+1;
        resnum = [resnum, atom.resSeq];
        chain  = [chain, atom.chainID];
        c5prime = [c5prime; atom.X, atom.Y, atom.Z];
        o5prime = [o5prime; NaN, NaN, NaN];
        c4prime = [c4prime; NaN, NaN, NaN];
    end
    if strcmp( atom.AtomName, 'O5''' )
        count = find( resnum == atom.resSeq );
        if ( ~isempty(count) ); o5prime(count,:) = [atom.X, atom.Y, atom.Z]; end
    end
    if strcmp( atom.AtomName, 'C4''' )
        count = find( resnum == atom.resSeq );
        if ( ~isempty(count) ); c4prime(count,:) = [atom.X, atom.Y, atom.Z]; end
    end
    xyz =[xyz; atom.X, atom.Y, atom.Z];
end


%%
%Set up coordinate frames
ctr = []; M = [];
for n = 1:size( o5prime, 1 );
    ctr(:,n) = o5prime(n,:)';
    M(:,:,n) = get_coordinate_frame( o5prime(n,:), c5prime(n,:),c4prime(n,:) );
end

%%
% draw it
cla
% coordinates only
%plot3( o5prime(:,1), o5prime(:,2),o5prime(:,3),'linew',2 ); hold on
%plot3( c5prime(:,1), c5prime(:,2),c5prime(:,3) )
%plot3( c4prime(:,1), c4prime(:,2),c4prime(:,3) )
%plot3( xyz(:,1),xyz(:,2),xyz(:,3),'o'); hold on
draw_trace( ctr, M );

%%
% get frame-to-frame transforms
t = []; % translations
R = []; % 3x3 rotation matrices
for n = 1:(size( o5prime, 1)-1);
    % later need to put in a filter for chainbreaks!
    [t(:,n),R(:,:,n)] = get_transform( ctr(:,n), M(:,:,n), ctr(:,n+1), M(:,:,n+1));
end


%% 
% Cool let's generate a random trajectory
N = 100;
[x,m] = get_random_trace(N, t, R );

%%
% Let's get some statistics, and estimate return probability (C_eff for
% circular RNA)

%which_N = [2:9 10:5:40];
N = 2; % number of steps

NITER = 1000;
[pts_x,pts_m] = get_pts_forward( NITER, N, t, R);

cla;plot3( pts_x(1,:),pts_x(2,:),pts_x(3,:),'o'); axis equal
% then collect histograms
     
  
% try KDE-based calc of C_eff
% convert 3x3 rotation matrices to angle-axis (Euler vector)
pts_EV=SpinCalc('DCMtoEV',pts_m,0,0);  % output is unit vector, angle in degrees
pts_EV3 = [];
% convert to axis vector v_x,v_y,v_z; with length equal to rotation angle in radians
for n = 1:size( pts_EV,1); pts_EV3(n,:) = pts_EV(n,1:3) * pts_EV(n,4) * pi/180.0; end;   

% points in 6D SE(3) space: x,y,z, v_x, v_y, v_z 
pts_f = [pts_x', pts_EV3];

s = get_kde_bandwidth( pts_f );
pts_r = [0,0,0,0,0,0];
p = mvksdensity(pts_f,pts_r,'Bandwidth',s)'

    
    
    
    
    
