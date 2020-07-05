function BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
% BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
%
% Get continguous dinucleotides from PDB.
%
% Input
%  pdbstruct   = PDB struct as read in from pdbread()
%
% Output
%  BB_dinucleotides = cell of structs with fields:
%                       resnum1,chain1,segid1, resnum2,chain2,segid2
%                     that correspond to contiguous dinucleotides.
%
% (C) R. Das, Stanford University, 2020

resnum = [];
chain = '';
segid = {};
o3prime = []; p = [];
count = 0;
for i = 1:length( pdbstruct.Model.Atom )
    atom = pdbstruct.Model.Atom(i);
    if strcmp( atom.AtomName, 'P' )
        count = count+1;
        resnum = [resnum, atom.resSeq];
        chain  = [chain, atom.chainID];
        segid{count}  = strip(atom.segID);
        p       = [p; atom.X, atom.Y, atom.Z];
        o3prime = [o3prime; NaN, NaN, NaN];
    end
    o3prime = fill_in_xyz(atom,resnum,chain,'O3''', o3prime );
end

% detect chainbreaks
BB_dinucleotides = {};
count = 0;
for n = 1:(size( p, 1 )-1);
    d = norm(o3prime(n,:) - p(n+1,:));
    if (d>2.0); continue; end
    count = count+1;
    y = struct();

    y.resnum1 = resnum(n);
    y.chain1 = chain(n);
    y.segid1 = segid{n};
    
    y.resnum2 = resnum(n+1);
    y.chain2 = chain(n+1);
    y.segid2 = segid{n+1};

    BB_dinucleotides = [BB_dinucleotides, y]; 
end

