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
% Cool let's generate a bunch of random trajectories
N = 100; % number of steps
x = zeros(3,N); m = zeros(3,3,N);
x(:,1) = [0,0,0]; % trajectory
m(:,:,1) = [1 0 0; 0 1 0; 0 0 1]; % orthonormal coordinate frame
ntransforms = size(t,2);
for n = 2:N
    j = randi(ntransforms);
    %j = n+20; %randi(ntransforms);
    x(:,n)  = x(:,n-1) + m(:,:,n-1)*t(:,j);
    m(:,:,n)= m(:,:,n-1)*R(:,:,j);
        
% Sanity check:
% ctr2 == ctr1 + M1*t

end
cla
draw_trace(x,m)



     
    
    
    
    
    
    
    
    
