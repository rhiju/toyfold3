function dinucleotides_reverse = reverse_dinucleotides( dinucleotides );
% dinucleotides_reverse = reverse_dinucleotides( dinucleotides );
% 
% Reverse resnum1/chain1/segid1 <--> resnum2/chain2/segid2 
%
% INPUT
%   dinucleotides =  cell of structs
%
% OUTPUT
%  dinucleotides_reverse = cell of structs with reversed 
%                           resnum1/chain1/segid1 <--> resnum2/chain2/segid2 
% (C) R. Das, Stanford, 2020

dinucleotides_reverse = {};

for i = 1:length( dinucleotides);
    x = dinucleotides{i};
    
    y = struct();
    y.resnum1 = x.resnum2;
    y.chain1 = x.chain2;
    y.segid1 = x.segid2;
    
    y.resnum2 = x.resnum1;
    y.chain2 = x.chain1;
    y.segid2 = x.segid1;
    
    dinucleotides_reverse = [dinucleotides_reverse, y];
end
