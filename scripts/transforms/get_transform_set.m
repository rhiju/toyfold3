function TransformSet = get_transform_set(pdbstruct,residue_pairs,start_triad,end_triad);
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
%  pdbstruct   = PDB struct as read in from pdbread()
%  residue_pairs = cell of structs with fields:
%                    resnum1,chain1,segid1, resnum2,chain2,segid2
%                  or string for type of transform ('BB')
%
%  start_triad = cell with 3 fields -- atom names to use for start triad
%  end_triad   = cell with 3 fields -- atom names to use for end triad
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
[ctr1,M1,resnum,chain,segid] = get_frames( pdbstruct, start_triad );
[ctr2,M2,resnum,chain,segid] = get_frames( pdbstruct, end_triad );
if ischar(residue_pairs) && strcmp(residue_pairs, 'BB' ); residue_pairs = get_BB_dinucleotides(pdbstruct); end
    
for n = 1:length( residue_pairs )
    respair = residue_pairs{n};
    i = intersect(find(resnum==respair.resnum1 & strcmp(segid, respair.segid1)), strfind(chain,respair.chain1));
    j = intersect(find(resnum==respair.resnum2 & strcmp(segid, respair.segid2)), strfind(chain,respair.chain2));
    if length(i)~= 1; continue; end;
    if length(j)~= 1; continue; end;
    count = count + 1;
    [t(:,count),R(:,:,count)] = get_transform( ctr1(:,i), M1(:,:,i), ctr2(:,j), M2(:,:,j));
end
TransformSet = struct( 't',t,'R',R);
TransformSet = fill_T6_from_t_and_R( TransformSet );


