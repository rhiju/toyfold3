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

%% Let's get an adenine riboswitch
pdbstruct_ade = pdbread( '../data/1y26.pdb');
kiss_dinucleotide = struct('resnum1',59,'chain1','X','segid1','','resnum2',40,'chain2','X','segid2','');
TransformLibrary.ade_kiss = get_transform_set( pdbstruct_ade, {kiss_dinucleotide},  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );


J23_dinucleotide = struct('resnum1',45,'chain1','X','segid1','','resnum2',54,'chain2','X','segid2','');
TransformLibrary.native_j23 = get_transform_set( pdbstruct_ade, {J23_dinucleotide},  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );

%% ((((((......[[.))))))........((((((]].....))))))
%                 |__________________|
%                40                  59
%
step_types = [{'ade_kiss'},repmat({'BB_stem1'},1,6-1),repmat({'BB'},1,8+1),repmat({'BB_stem1'},1,6-1)];
% circshift to make sampling easier -- final close in loop.
step_types = [repmat({'BB'},1,4),repmat({'BB_stem1'},1,6-1),{'ade_kiss'},repmat({'BB_stem1'},1,6-1),repmat({'BB'},1,5)]; 
sample_motif( step_types, 1000, TransformLibrary );

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling the 3WJ 'loop'
% first as just a free 8-mer.
NITER = 1000; clear loop8_traces;
for i = 1:NITER; loop8_traces{i} = get_random_trace(repmat({'BB'},1,8+1), TransformLibrary, 0); end
TransformLibrary.loop8mer = get_transforms_from_traces(loop8_traces,1,1+8+1);

%% check sampling with this library
step_types = [repmat({'BB_stem'},1,6-1),{'ade_kiss'},repmat({'BB_stem'},1,6-1),'loop8mer']; 
%step_types = [{'ade_kiss'},repmat({'BB_stem1'},1,6-1),'loop8mer',repmat({'BB_stem1'},1,6-1)]; 
sample_motif( step_types, 1000, TransformLibrary );

%% Now sample loop within context of 3WJ
step_types = [repmat({'BB'},1,8+1),'BP',repmat({'BB'},1,2+1),'BP',repmat({'BB'},1,3+1),'BP'];
%step_types = [repmat({'BB'},1,1+1),'BP',repmat({'BB'},1,1+1),'BP',repmat({'BB'},1,1+1),'BP']
loop8_3WJ_traces = sample_circle_trajectories(step_types, TransformLibrary, 500);
%%
for i = 1:4; subplot(2,2,i); draw_trace( loop8_3WJ_traces{i}, step_types ); end;
TransformLibrary.loop8mer_3WJ = get_transforms_from_traces(loop8_3WJ_traces,1,1+8+1);

%% test length matching.
%
%   A         B
% xxxxx     xxxxxx
% |||||  8  ||||||
% xxxxxxxxxxxxxxxx
% |______________|
%    C_eff
%
% Fix total length at L.
%  Scan A and B
%  Expect  a 'resonance' when A = B = 6, which corresponds
%  to 1y26 adenine riboswitch.
figure(1);
NITER = 500;
which_lengths = [1:30];
C_eff_matrix = [];
BB_stem_type = 'BB_stem';
%step_types_loop = repmat({'BB'},1,8+1); step_types_tert = 'ade_kiss';
%step_types_loop = 'loop8mer'; step_types_tert = 'ade_kiss';
%step_types_loop = 'loop8mer_3WJ'; step_types_tert = 'ade_kiss';
step_types_loop = 'native_j23'; step_types_tert = 'ade_kiss';
for i = 1:length(which_lengths)
    for j = 1:length(which_lengths)
        if (i <= size(C_eff_matrix,1) &  j <= size(C_eff_matrix,2) & C_eff_matrix(i,j)>0)  continue; end;
        A = which_lengths(i);
        B = which_lengths(j);
        fprintf( 'Doing %d,%d...\n',A,B);
        %step_types = [{'ade_kiss'},repmat({BB_stem_type},1,A-1),step_types_loop,repmat({BB_stem_type},1,B-1)];
        step_types = [repmat({BB_stem_type},1,A-1),step_types_loop,repmat({BB_stem_type},1,B-1),{step_types_tert}];
        [C_eff_matrix(i,j),C_eff_matrix_err(j,j)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
    end
end

%%
clf
imagesc(which_lengths,which_lengths,C_eff_matrix)
xlabel( 'length of helix A'); 
ylabel( 'length of helix B'); 
h=title( sprintf('A-8mer-B, from L2-L3 kissing loop (ade). BB stem type = %s, loop type = %s',BB_stem_type,step_types_loop ) )
set(h,'interp','none')
set(gcf, 'PaperPositionMode','auto','color','white');
set(gca,'ydir','normal');
colorbar;
rectangle('position',[5.5 5.5 1 1]);

