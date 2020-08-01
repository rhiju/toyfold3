%%
pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
stems = read_stems_toyfold3( '../data/4ybb_DIII.pdb.stems.txt' );
%pdbstruct = pdbread( '../data/4ybb_23S.pdb');
%stems = read_stems_toyfold3( '../data/4ybb_23S.pdb.stems.txt' );
%%
tic
BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
TransformLibrary.BB = get_transform_set( pdbstruct, BB_dinucleotides, {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc
%%
tic
base_pairs = get_base_pairs_from_stems_toyfold3( stems );
TransformLibrary.BP = get_transform_set( pdbstruct, base_pairs,  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc

%%
pdbstruct = pdbread( '../data/2gjw_nuclease_inline_conf.pdb');
TransformLibrary.Inline = get_transform_set( pdbstruct, 'BB', {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's take a look at closed trajectories
[C_eff_BPS,C_eff_BPS] = sample_motif( {'BP','BB','BP','BB'}, 1000, TransformLibrary );

%%
[C_eff_singlentbulge,C_eff_singlentbulge_err] = sample_motif( {'BB','BP','BB','BP','BB'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge,C_eff_singlentbulge_x_err] = sample_motif( {'BP','BB','BP','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge_inline0,C_eff_singlentbulge_inline0_err] = sample_motif( {'BB','BP','BB','BP','Inline'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge_inline0x,C_eff_singlentbulge_inline0x_err] = sample_motif( {'BP','BB','BP','Inline','BB'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge_inline1,C_eff_singlentbulge_inline1_err] = sample_motif( {'Inline','BP','BB','BP','BB'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge_inline1x,C_eff_singlentbulge_inline1x_err] = sample_motif( {'BP','BB','BP','BB','Inline'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge_inline1xx,C_eff_singlentbulge_inline1xx_err] = sample_motif( {'BP','BB','Inline','BP','BB',}, 1000, TransformLibrary );


%%
[C_eff_tetraloop,C_eff_tetraloop_err] = sample_motif( {'BP','BB','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_tetraloop_inline0,C_eff_tetraloop_inline0_err] = sample_motif( {'BP','Inline','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_tetraloop_inline1,C_eff_tetraloop_inline1_err] = sample_motif( {'BP','BB','Inline','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_tetraloop_inline2,C_eff_tetraloop_inline2_err] = sample_motif( {'BP','BB','BB','Inline','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_tetraloop_inline3,C_eff_tetraloop_inline3_err] = sample_motif( {'BP','BB','BB','BB','Inline','BB'}, 1000, TransformLibrary );
%%
[C_eff_tetraloop_inline4,C_eff_tetraloop_inline4_err] = sample_motif( {'BP','BB','BB','BB','BB','Inline'}, 1000, TransformLibrary );
%%


%%
figure(1)
[C_eff_triloop,C_eff_triloop_err] = sample_motif( {'BP','BB','BB','BB','BB'}, 1000, TransformLibrary );
figure(2)
[C_eff_triloop_inline0,C_eff_triloop_inline0_err] = sample_motif( {'BP','Inline','BB','BB','BB'}, 1000, TransformLibrary );
figure(3)
[C_eff_triloop_inline1,C_eff_triloop_inline1_err] = sample_motif( {'BP','BB','Inline','BB','BB'}, 1000, TransformLibrary );
figure(4)
[C_eff_triloop_inline2,C_eff_triloop_inline2_err] = sample_motif( {'BP','BB','BB','Inline','BB'}, 1000, TransformLibrary );
figure(5)
[C_eff_triloop_inline3,C_eff_triloop_inline3_err] = sample_motif( {'BP','BB','BB','BB','Inline'}, 1000, TransformLibrary );
%%


%% long loop (decaloop)
figure(1)
[C_eff_decaloop] = sample_motif( {'BP','BB','BB','BB','BB','BB','BB','BB','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_decaloop_inline5] = sample_motif( {'BP','BB','BB','BB','BB','BB','Inline','BB','BB','BB','BB','BB'}, 1000, TransformLibrary );


% Does 3WJ accommodate inline? Hammerhead mimic
%%
[C_eff_3WJ_tetraloop] = sample_motif( {'BB','BP','BB','BP','BB','BP','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_3WJ_tetraloop_inlineA] = sample_motif( {'BB','BP','Inline','BP','BB','BP','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_3WJ_tetraloop_inlineB] = sample_motif( {'BB','BP','BB','BP','Inline','BP','BB','BB','BB','BB'}, 1000, TransformLibrary );


% Try to really mimic hammerhead (3zp8):
%%
[C_eff_HH] = sample_motif( {'BB','BP','BB','BB','BP','BB','BP','BB','BB','BP','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_HH_inline] = sample_motif( {'BB','BP','BB','Inline','BP','BB','BP','BB','BB','BP','BB','BB','BB'}, 1000, TransformLibrary );



% How about a 4WJ?
%%
[C_eff_4WJ] = sample_motif( {'BB','BP','BB','BP','BB','BP','BB','BP'}, 1000, TransformLibrary );
%%
[C_eff_4WJ_Inline] = sample_motif( {'BB','BP','Inline','BP','BB','BP','BB','BP'}, 1000, TransformLibrary );

