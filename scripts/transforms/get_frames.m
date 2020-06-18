function [ctr,M,ok] = get_frames( pdbstruct, atom_triad );
%
% Grab coordinate frames at O5'
%
%  Would be better to go through HETATMS too and
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
%  atom_triad = (cell of strings) Triad of atom names, and extra flag
%                    {'C5''','C4''','C3'''} or
%                    {'C5''','C4''','C3''',1} to get to 'next base', and check if no chainbreak. 
%
% Outputs
%  ctr = [3 x N] coordinates of the trace
%  M   = [3 x 3 x N] orthonormal frames of the trace
%  ok  = flag that outputs special information, e.g., continuous chain. 
%
% (C) R. Das, Stanford 2020

count = 0;
resnum = [];
chain = '';
atom_name1 = atom_triad{1};
atom_name2 = atom_triad{2};
atom_name3 = atom_triad{3};
next_base = 0;
if length( atom_triad ) == 4 && atom_triad{4} == 1; next_base = 1; end;

atom1 = []; atom2 = []; atom3 = []; 
o3prime = []; p = [];
for i = 1:length( pdbstruct.Model.Atom )
    atom = pdbstruct.Model.Atom(i);
    if strcmp( atom.AtomName, 'P' )
        count = count+1;
        resnum = [resnum, atom.resSeq];
        chain  = [chain, atom.chainID];
        p       = [p; atom.X, atom.Y, atom.Z];
        o3prime = [o3prime; NaN, NaN, NaN];
        atom1 = [atom1; NaN, NaN, NaN];
        atom2 = [atom2; NaN, NaN, NaN];
        atom3 = [atom3; NaN, NaN, NaN];
    end
    atom1 = fillin(atom,resnum,chain,atom_name1, atom1 );
    atom2 = fillin(atom,resnum,chain,atom_name2, atom2 );
    atom3 = fillin(atom,resnum,chain,atom_name3, atom3 );
    o3prime = fillin(atom,resnum,chain,'O3''', o3prime );
end

%%
%Set up coordinate frames
ctr = []; M = [];
for n = 1:size( atom1, 1 );
    ctr(:,n) = atom1(n,:)';
    M(:,:,n) = get_coordinate_frame( atom1(n,:), atom2(n,:), atom3(n,:) );
end

ok = ones(1,size(atom1,1));

if next_base
    % detect chainbreaks
    chainbreak = zeros( 1, size( atom1, 1 ) );
    chainbreak(end) = 1;
    for n = 1:(size( atom1, 1 )-1);
        d = norm(o3prime(n,:) - p(n+1,:));
        if (d>2.0); chainbreak(n) = 1;  end
    end

    ctr(:,1:end-1) = ctr(:,2:end);
    M(:,:,1:end-1) = M(:,:,2:end);
    ctr(:,end) = NaN;
    M(:,:,end) = NaN;
    ok = 1-chainbreak;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function xyz = fillin(atom,resnum,chain,atomname, xyz );
if strcmp( atom.AtomName, atomname )
    count = find( resnum == atom.resSeq & strfind(chain,atom.chainID) );
    if ( ~isempty(count) ); 
        xyz(count,:) = [atom.X, atom.Y, atom.Z]; 
    end
end
