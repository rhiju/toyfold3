function  BB_dinucleotides = get_BB_from_stems( stems )
% BB_dinucleotides = get_BB_from_stems_toyfold3( stems )
%
%  Infer (i,i+1) from each side of stems
% 
% (C) R. Das, Stanford University, 2020
BB_dinucleotides = {};
for n = 1:length( stems )
    stem = stems{n};
    N = length( stem.resnum1 );
    for k = 1:(N-1)
        y = struct();
        y.resnum1 = stem.resnum1(k);
        y.chain1  = stem.chain1(k);
        y.segid1  = stem.segid1{k};
        y.resnum2 = stem.resnum1(k+1);
        y.chain2  = stem.chain1(k+1);
        y.segid2  = stem.segid1{k+1};
        BB_dinucleotides = [ BB_dinucleotides, {y} ];
    end
    for k = 1:(N-1)
        y = struct();
%         y.resnum1 = stem.resnum2(N-k+1);
%         y.chain1  = stem.chain2(N-k+1);
%         y.segid1  = stem.segid2{N-k+1};
%         y.resnum2 = stem.resnum2(N-k);
%         y.chain2  = stem.chain2(N-k);
%         y.segid2  = stem.segid2{N-k};
        y.resnum1 = stem.resnum2(k);
        y.chain1  = stem.chain2(k);
        y.segid1  = stem.segid2{k};
        y.resnum2 = stem.resnum2(k+1);
        y.chain2  = stem.chain2(k+1);
        y.segid2  = stem.segid2{k+1};
        BB_dinucleotides = [ BB_dinucleotides, {y} ];
    end
end
