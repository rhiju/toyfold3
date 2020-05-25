function [ctr,M,chainbreak] = get_frames( pdbstruct );
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
% Outputs
%  ctr = [3 x N] coordinates of the trace
%  M   = [3 x 3 x N] orthonormal frames of the trace
%
% (C) R. Das, Stanford 2020

count = 0;
resnum = [];
chain = '';
o5prime = []; c5prime = []; c4prime = []; c3prime = []; o3prime = []; p = [];
for i = 1:length( pdbstruct.Model.Atom )
    atom = pdbstruct.Model.Atom(i);
    if strcmp( atom.AtomName, 'P' )
        count = count+1;
        resnum = [resnum, atom.resSeq];
        chain  = [chain, atom.chainID];
        p       = [p; atom.X, atom.Y, atom.Z];
        c5prime = [c5prime; NaN, NaN, NaN];
        o5prime = [o5prime; NaN, NaN, NaN];
        c4prime = [c4prime; NaN, NaN, NaN];
        c3prime = [c3prime; NaN, NaN, NaN];
        o3prime = [o3prime; NaN, NaN, NaN];
    end
    o5prime = fillin(atom,resnum,chain,'O5''', o5prime );
    c5prime = fillin(atom,resnum,chain,'C5''', c5prime );
    c4prime = fillin(atom,resnum,chain,'C4''', c4prime );
    c3prime = fillin(atom,resnum,chain,'C3''', c3prime );
    o3prime = fillin(atom,resnum,chain,'O3''', o3prime );
end

%%
%Set up coordinate frames
ctr = []; M = [];
for n = 1:size( o5prime, 1 );
    ctr(:,n) = c5prime(n,:)';
    M(:,:,n) = get_coordinate_frame( c5prime(n,:), c4prime(n,:), c3prime(n,:) );
end

% detect chainbreaks
chainbreak = zeros( 1, size( o5prime, 1 ) );
chainbreak(end) = 1;
for n = 1:(size( o5prime, 1 )-1);
    d = norm(o3prime(n,:) - p(n+1,:));
    if (d>2.0); chainbreak(n) = 1;  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function xyz = fillin(atom,resnum,chain,atomname, xyz );
if strcmp( atom.AtomName, atomname )
    count = find( resnum == atom.resSeq & strfind(chain,atom.chainID) );
    if ( ~isempty(count) ); 
        xyz(count,:) = [atom.X, atom.Y, atom.Z]; 
    end
end
