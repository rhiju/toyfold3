function motifs = get_unique_motifs(Nloop, Lmax, motifs);
% motifs = get_unique_motifs(Nloop, Lmax);
%
% Enumerate motifs (e.g., 3-way junctions).
%
% Inputs
%  Nloop = number of loops (or stems). E.g., 3 for three-way-junction.
%  Lmax  = maximum number of loop nucleotides, summing loops all the
%           way around junction.
%  motifs = starting cell of motifs [Optional, default is empty]
%
% Output
%  motifs = cell of unique motifs represented as arrays; each array
%            gives loop lengths.
%
% (C) R. Das, Stanford University, 2020

if ~exist( 'motifs', 'var') motifs = {}; end;
for L = 0:Lmax
    for r = 1:(Nloop-1)^(L+1) % total number of options for lengths
        new_motif = [];
        count = r-1;
        for i = 1:(Nloop-1)
            new_motif = [new_motif, mod(count,L+1) ];
            count = floor(count/(L+1));
        end
        l = L-sum(new_motif);
        if ( l < 0 ); continue; end;
        new_motif = [new_motif, l];
        
        ok = 1;
        for q = 1:length(motifs)
            if length(motifs{q})~=Nloop; continue;end
            for n = 1:Nloop
                if all(circshift(new_motif,n)==motifs{q}); ok = 0; break; end;
            end
        end
        if ~ok; continue;end;
        
        motifs = [ motifs, {new_motif}];
    end
end