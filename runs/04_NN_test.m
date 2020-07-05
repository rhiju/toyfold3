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
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cool let's generate a random trajectory
%% base pair step, a.k.a. stacked pair  [no closure]
step_types = {'BB','BP','BB','BP'}; 
x = get_random_trace(step_types, TransformLibrary );
%% 1x1 internal loop [no closure]
step_types = {'BB','BB','BP','BB','BB','BP'}; 
x = get_random_trace(step_types, TransformLibrary );


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's take a look at closed trajectories
[C_eff_BPS,C_eff_BPS] = sample_motif( {'BP','BB','BP','BB'}, 1000, TransformLibrary );
%%
[C_eff_singlentbulge,C_eff_singlentbulge_err] = sample_motif( {'BB','BP','BB','BP','BB'}, 1000, TransformLibrary );
%%
[C_eff_tetraloop,C_eff_tetraloop] = sample_motif( {'BP','BB','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_3WJ,C_eff_3WJ] = sample_motif( {'BB','BP','BB','BP','BB','BP'}, 1000, TransformLibrary );
%%
[C_eff_3WJ_tetraloop,C_eff_3WJ_tetraloop] = sample_motif( {'BB','BP','BB','BP','BB','BP','BB','BB','BB','BB'}, 1000, TransformLibrary );
%%
[C_eff_4WJ,C_eff_4WJ] = sample_motif( {'BB','BP','BB','BP','BB','BP','BB','BP'}, 1000, TransformLibrary );

%% Scan through circular BB lengths,hairpin lengths,
loop_lengths_BB = [2:20]; NITER = 5000;
out_BB = scan_loop_length(loop_lengths_BB,NITER,TransformLibrary);

%%
loop_lengths = [1:20];
NITER = 2000;
out_HP = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB'});
out_0xN = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB','BP','BB'});
out_1xN = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB','BB','BP','BB'});
out_2xN = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB','BB','BB','BP','BB'});
out_0x0xN = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB','BP','BB','BP','BB'});
out_0x0x0xN = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB','BP','BB','BP','BB','BP','BB'});

%%
clf
errorbar( 0, C_eff_overlap_f_BPS(2), C_eff_overlap_f_BPS_err(2),'linew',2); hold on
errorbar( loop_lengths, out_0xN.C_eff,  out_0xN.C_eff_err,'linew',2); hold on
errorbar( loop_lengths, out_1xN.C_eff,  out_1xN.C_eff_err,'linew',2); hold on
errorbar( loop_lengths, out_2xN.C_eff,  out_2xN.C_eff_err,'linew',2); hold on
errorbar( loop_lengths, out_0x0x0xN.C_eff,  out_0x0x0xN.C_eff_err,'linew',2); hold on
errorbar( loop_lengths, out_HP.C_eff,  out_HP.C_eff_err,'linew',2); hold on
errorbar( loop_lengths, out_0x0xN.C_eff,  out_0x0xN.C_eff_err,'linew',2); hold on
errorbar( loop_lengths_BB, out_BB.C_eff,  out_BB.C_eff_err,'linew',2); hold on
set(gca,'fontweight','bold','yscale','log');xlabel('N');ylabel('C_{eff} (M)');
legend( 'Base pair stack','0xN internal (asymmetric bulge)','1xN internal','2xN internal','0x0x0xN (4WJ)',...
    'Hairpin','0x0xN (3WJ)','Circular BB');
set(gcf, 'PaperPositionMode','auto','color','white');
ylim([1e-6 1e3]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kind of advanced --> Let's sample hairpins and check that answer
%%  is independent of which link has the inline conformation and
%%  also how many links we go through to get overlap prob.
%N = 4; NITER = 2000;
N = 10; NITER = 1000;
clear all_pts_f all_pts_r C_eff_overlap_f C_eff_overlap_r;
cla;
legends = {};
% hairpins -- check independence w.r.t. BP location within 'motif'.
for j = 1:(N+2);
    step_types = repmat({'BB'},1,N+2);
    step_types{j} = 'BP';
    step_type_sets{j} = step_types;
    [all_pts_f{j}, all_pts_r{j}] = get_all_pts( step_types, NITER, TransformLibrary );
    [C_eff_overlap_f{j}, C_eff_overlap_r{j}] = get_C_eff_overlap( all_pts_f{j}, all_pts_r{j} );
end
%%
clf
for j = 1:(N+2);
    plot( 1:(N+2),C_eff_overlap_f{j} ); hold on
    plot( 1:(N+2),C_eff_overlap_r{j}); hold on
    set(gca,'fontweight','bold','xgrid','on'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
    title(sprintf('Forward/reverse overlap molarity for hairpin with loop length %d',N) );
    legends = [legends, {sprintf('C_eff_overlap_f, BP link at %d',j),sprintf('C_eff_overlap_r, BP link at %d',j)}];
end
h = legend( legends); set(h,'interpreter','none');

%% Let's see if circularized trajectories make sense
j = N/2;
sample_circle_trajectory( step_type_sets{j}, TransformLibrary, all_pts_r{j} );


