function TransformSet = get_transform_set(pdbstruct,start_triad,end_triad);
%
% get frame-to-frame transforms
%
% Grab coordinate frames at O5'
%
%  Would be better to go through HETATMS too and
%   figure out order in some sane way?
%
%    O5' <-- center
%     z\      z points from O5' to C5'
%       C5'   x points perpendicular to z towards C4'
%       |    [y completes orthonormal frame.]
%       C4'
%
% Inputs
%  pdbstruct = PDB struct as read in from pdbread()
%
% Output
%  TransformSet = struct with two fields:
%     t = [3 x Nframes] library of translations from nt to nt. 
%     R = [3 x 3 x Nframes] library of rotations from nt to nt. 
%
% (C) R. Das, Stanford 2020

t = []; % translations
R = []; % 3x3 rotation matrices
count = 0;
[ctr1,M1,ok1] = get_frames( pdbstruct, start_triad );
[ctr2,M2,ok2] = get_frames( pdbstruct, end_triad );

for n = find( ok1 & ok2)
    % later need to put in a filter for chainbreaks!
    count = count+1;
    [t(:,count),R(:,:,count)] = get_transform( ctr1(:,n), M1(:,:,n), ctr2(:,n), M2(:,:,n));
end
TransformSet = struct( 't',t,'R',R);

