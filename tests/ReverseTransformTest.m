pdbstruct = pdbread( 'example_data/4ybb_DIII_just3nts.pdb');
start_triad = {'C5''','C4''','C3'''};
[ctr,M] = get_frames( pdbstruct, start_triad );
ctr1 = ctr(:,1);
ctr2 = ctr(:,2);
M1 = M(:,:,1);
M2 = M(:,:,2);

[tf,Rf] = get_transform( ctr1, M1, ctr2, M2);
[tr,Rr] = get_transform( ctr2, M2, ctr1, M1);
[tf_rev, Rf_rev] =  reverse_transform( tf,Rf);

assert( all(abs(tf_rev - tr)<1e-5) )
assert( all(abs(Rf_rev(:) - Rr(:))<1e-5) )