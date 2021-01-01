%%
pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
stems = read_stems_toyfold3( '../data/4ybb_DIII.pdb.stems.txt' );
%pdbstruct = pdbread( '../data/4ybb_23S.pdb');
%stems = read_stems_toyfold3( '../data/4ybb_23S.pdb.stems.txt' );
%%
tic
TransformLibary = struct();
BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
TransformLibrary.BB = get_transform_set( pdbstruct, BB_dinucleotides, {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc
%%
tic
base_pairs = get_base_pairs_from_stems_toyfold3( stems );
TransformLibrary.BP = get_transform_set( pdbstruct, base_pairs,  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );

%%  
% New ... needed to construct helices.
BB_stem_dinucleotides = get_BB_from_stems( stems );
TransformLibrary.BB_stem = get_transform_set( pdbstruct, BB_stem_dinucleotides,  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc

%% Just pick a single representative base pair step (since some of the dinucleotides above are a bit wacky)
TransformLibrary.BB_stem1 = get_transform_set( pdbstruct, BB_stem_dinucleotides(1),  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualize transforms -- should be pretty regular
%% C_eff matrix for
%     A         B
%   xxxxx     xxxx
% p |||||  q  ||||r
% xxxxxxxxxxxxxxxxx
% |_______________|
%    C_eff
%
% So total length = p+A+q+B+r
%

step_types = {};
p = 2; q = 4; r = 2;
A = 22; B = 22;
NITER = 1000;
step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r)];
[C_eff(1),C_eff_err(1)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r)];
[C_eff(2),C_eff_err(2)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r)];
[C_eff(3),C_eff_err(3)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r)];
[C_eff(4),C_eff_err(4)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );

fprintf( '\n');
tags = {'--','A-','-B','AB'};
for i = 1:4
    fprintf('%s %8.5f\n',tags{i},C_eff(i))
end
%%
NITER = 500;
%p = 2; q = 4; r = 2; which_lengths = [1:50];
p = 1; q = 1; r = 1; which_lengths = [1:50];
C_eff_scan = [];
for k = 1:length(which_lengths)
    A = which_lengths(k);
    B = A;    
    fprintf( 'Doing length %d...\n',which_lengths(k));
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r)];
    [C_eff_scan(k,1),C_eff_scan_err(k,1)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r)];
    [C_eff_scan(k,2),C_eff_scan_err(k,2)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r)];
    [C_eff_scan(k,3),C_eff_scan_err(k,3)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r)];
    [C_eff_scan(k,4),C_eff_scan_err(k,4)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
end

semilogy( which_lengths, C_eff_scan,'linew',2 );
legend(tags)
xlabel( 'length of helix A/B.'); 
ylabel( 'C_{eff} (M)')
title( sprintf('Vary helix length: p (5'' linker) = %d, q (helix-to-helix linker)= %d, r (3'' linker) = %d',p,q,r ) )
set(gcf, 'PaperPositionMode','auto','color','white');

%%

