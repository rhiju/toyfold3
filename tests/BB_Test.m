pdbstruct = pdbread( 'example_data/4ybb_DIII_just3nts.pdb');
BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
assert( length(BB_dinucleotides) == 2);

pdbstruct = pdbread( 'example_data/4ybb_DIII_just3nts_skipmiddle.pdb');
BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
assert( length(BB_dinucleotides) == 0);
