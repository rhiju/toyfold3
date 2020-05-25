function [ctr,M] = get_frames( pdbstruct );
%
% Grab coordinate frames at O5'
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
% Outputs
%  ctr = [3 x N] coordinates of the trace
%  M   = [3 x 3 x N] orthonormal frames of the trace
%
% (C) R. Das, Stanford 2020

count = 0;
resnum = [];
chain = '';
o5prime = []; c5prime = []; c4prime = [];
xyz = []; % collect all xyz, just for checks
for i = 1:length( pdbstruct.Model.Atom )
    atom = pdbstruct.Model.Atom(i);
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