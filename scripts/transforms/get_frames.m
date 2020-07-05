function [ctr,M,resnum,chain,segid] = get_frames( pdbstruct, atom_triad );
%
% Grab coordinate frames at specified triad
%
% Assumes P is defined for each atom -- would be better to not assume that!
%
% Would be better to go through HETATMS too and
%   figure out order in some sane way?
%
%    C5' <-- center
%     z\      z points from C5' to C4'
%       C4'   x points perpendicular to z towards C3'
%       |    [y completes orthonormal frame.]
%       C3'
%
% Inputs
%  pdbstruct = PDB struct as read in from pdbread()
%  atom_triad = (cell of strings) Triad of atom names, e.g.
%                    {'C5''','C4''','C3'''}
%
% Outputs
%  ctr = [3 x N] coordinates of the trace
%  M   = [3 x 3 x N] orthonormal frames of the trace
%  resnum = [number array of length N] residue number for triad
%  chain  = [string of length N] chain for triad 
%  segid  = [cell of strings of length N] segID for triad
%
% (C) R. Das, Stanford 2020

count = 0;
resnum = [];
chain = '';
segid = {};
atom_name1 = atom_triad{1};
atom_name2 = atom_triad{2};
atom_name3 = atom_triad{3};

p = [];
atom1 = []; atom2 = []; atom3 = []; 
for i = 1:length( pdbstruct.Model.Atom )
    atom = pdbstruct.Model.Atom(i);
    if strcmp( atom.AtomName, 'P' )
        count = count+1;
        resnum = [resnum, atom.resSeq];
        chain  = [chain, atom.chainID];
        segid  = [segid, {strip(atom.segID)}];
        p       = [p; atom.X, atom.Y, atom.Z];
        atom1 = [atom1; NaN, NaN, NaN];
        atom2 = [atom2; NaN, NaN, NaN];
        atom3 = [atom3; NaN, NaN, NaN];
    end
    atom1 = fill_in_xyz(atom,resnum,chain,atom_name1, atom1 );
    atom2 = fill_in_xyz(atom,resnum,chain,atom_name2, atom2 );
    atom3 = fill_in_xyz(atom,resnum,chain,atom_name3, atom3 );
end

%%
%Set up coordinate frames
ctr = []; M = [];
for n = 1:size( atom1, 1 );
    ctr(:,n) = atom1(n,:)';
    M(:,:,n) = get_coordinate_frame( atom1(n,:), atom2(n,:), atom3(n,:) );
end
