function xyz = fill_in_xyz(atom,resnum,chain,atomname, xyz );
% xyz = fill_in_xyz(atom,resnum,chain,atomname, xyz );
%
% Helper function to fill in xyz based on resnum, chain, and atomname
%    discovered so far.
%
% (C) R. Das, Stanford University
if strcmp( atom.AtomName, atomname )
    count = find( resnum == atom.resSeq & strfind(chain,atom.chainID) );
    if ( ~isempty(count) ); 
        xyz(count,:) = [atom.X, atom.Y, atom.Z]; 
    end
end
