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
% Scan length of helix -- any regime with strong XOR?
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
% try an AND gate (based on length matching)
%
%     A         B      C
%   xxxx     xxxx  xxxxxxxx
% p ||||  q  ||||r |||||||| s
% xxxxxxxxxxxxxxxxxxxxxxxxxxx
% |_______________________|
%    C_eff
%
%  IF C = 20 bp, hope that A = B = 10 bp might
%  lead to AND behavior?
%

NITER = 500;
%p = 2; q = 4; r = 2; which_lengths = [1:50];
p = 1; q = 0; r = 1; s = 1; which_lengths = [5:5:50];
C_eff_scan_AND = [];
for k = 1:length(which_lengths)
    A = which_lengths(k);
    B = A; C = 2*A;    
    fprintf( 'Doing length %d...\n',which_lengths(k));
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r+1),repmat({'BB_stem'},1,C-1),repmat({'BB'},1,s)]; 
    [C_eff_scan_AND(k,1),C_eff_scan_AND_err(k,1)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r+1),repmat({'BB_stem'},1,C-1),repmat({'BB'},1,s)];
    [C_eff_scan_AND(k,2),C_eff_scan_AND_err(k,2)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r+1),repmat({'BB_stem'},1,C-1),repmat({'BB'},1,s)];
    [C_eff_scan_AND(k,3),C_eff_scan_AND_err(k,3)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r+1),repmat({'BB_stem'},1,C-1),repmat({'BB'},1,s)];
    [C_eff_scan_AND(k,4),C_eff_scan_AND_err(k,4)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
end

semilogy( which_lengths, C_eff_scan_AND,'linew',2 );
legend(tags)
xlabel( 'length of helix A/B.'); 
ylabel( 'C_{eff} (M)')
title( sprintf('AND -- Vary helix length: p (5'' linker) = %d, q (helix-to-helix linker)= %d, r (helix-to-const-helix) = %d, s (3'' linker) = %d',p,q,r,s ) )
set(gcf, 'PaperPositionMode','auto','color','white');


%%
% try to test length matching.
%
%     A         B
%   xxxxx     xxxx
% p ||||| l q ||||r
% xxxxxxxxxxxxxxxxx
% |_______________|
%    C_eff
%
% So total length = p+A+q+B+r
%  Fix B's length.
%  Scan A, and set l = B - A (to keep total loop length constant). 
%  Expect  a 'resonance' when A = B.


NITER = 500;
%B = 20; p = 1; q = 0; r = 1; which_lengths = [1:B];
%B = 20; p = 2; q = 4; r = 2; which_lengths = [1:B];
B = 50; p = 2; q = 4; r = 2; which_lengths = [1:B];

C_eff_scanA = zeros(length(which_lengths),4);
for k = 1:length(which_lengths)
    A = which_lengths(k);
    l = B - A;
    fprintf( 'Doing length %d...\n',which_lengths(k));
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,l+q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r)];
    [C_eff_scanA(k,1),C_eff_scanA_err(k,1)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,l+q+1),repmat({'BB'     },1,B-1),repmat({'BB'},1,r)];
    [C_eff_scanA(k,2),C_eff_scanA_err(k,2)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB'     },1,A-1),repmat({'BB'},1,l+q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r)];
    [C_eff_scanA(k,3),C_eff_scanA_err(k,3)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'},1,A-1),repmat({'BB'},1,l+q+1),repmat({'BB_stem'},1,B-1),repmat({'BB'},1,r)];
    [C_eff_scanA(k,4),C_eff_scanA_err(k,4)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
end


semilogy( which_lengths, C_eff_scanA,'linew',2 );
legend(tags)
xlabel( 'length of helix A'); 
ylabel( 'C_{eff} (M)')
title( sprintf('p-A-l-q-B-r.  p = %d, A = varies, l=B-A (keep loop length fixed), q = %d, B = %d, r  = %d',p,q,B,r ) )
set(gcf, 'PaperPositionMode','auto','color','white');



%% try to test length matching.
%
%     A         B
%   xxxxx     xxxx
% p ||||| l   ||||r
% xxxxxxxxxxxxxxxxx
% |_______________|
%    C_eff
%
% Fix total length at L.
%  Scan A and B, set l so that p+A+l+B+r = L
%  Expect  a 'resonance' when A = B.

NITER = 500;
%L = 12; p = 1; r = 1; which_lengths = [0:1:10];
%L = 42; p = 1; r = 1; which_lengths = [0:1:40];
L = 102; p = 1; r = 1; which_lengths = [0:1:100];
%C_eff_matrix = [];
for i = 1:length(which_lengths)
    for j = 1:length(which_lengths)
        if (i <= size(C_eff_matrix,1) &  j <= size(C_eff_matrix,2) & C_eff_matrix(i,j)>0)  continue; end;
        A = which_lengths(i);
        B = which_lengths(j);
        l = L-(A+B+p+r);
        if l<0; continue; end
        fprintf( 'Doing %d,%d...\n',A,B);
        step_types = [{'BP'},repmat({'BB'},1,p),repmat({'BB_stem'     },1,A-1),repmat({'BB'},1,l+1),repmat({'BB_stem'     },1,B-1),repmat({'BB'},1,r)];
        [C_eff_matrix(i,j),C_eff_matrix_err(j,j)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    end
end

imagesc(which_lengths,which_lengths,C_eff_matrix)
xlabel( 'length of helix A'); 
ylabel( 'length of helix B'); 
title( sprintf('p-A-l-B-r. Total L=%d fixed, p = %d,  r  = %d',L,p,r ) )
set(gcf, 'PaperPositionMode','auto','color','white');
set(gca,'ydir','normal');
colorbar;