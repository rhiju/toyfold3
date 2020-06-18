%% Let's get the inline conformation into our library:
load toyfold3_test.mat TransformLibrary
pdbstruct = pdbread('4ybb_23S.pdb');
TransformSet = get_transform_set( pdbstruct, {'C5''','C4''','C3'''},{'C5''','C4''','C3''',1} );
TransformLibrary.BB = TransformSet;

pdbstruct = pdbread( '../data/2gjw_nuclease_inline_conf.pdb');
TransformSet = get_transform_set( pdbstruct, {'C5''','C4''','C3'''},{'C5''','C4''','C3''',1} );
TransformLibrary.Inline = TransformSet;

%% Let's sample trajectories (not circles) and check that answer
%%  is independent of which link has the inline conformation and
%%  also how many links we go through to get overlap prob.
NITER = 1000;
N = 10; % N = 4;
clear all_pts_f all_pts_r C_eff_overlap_f C_eff_overlap_r;
cla;
legends = {};
for j = 1:N;
    step_types = repmat({'BB'},1,N);
    step_types{j} = 'Inline';
    step_type_sets{j} = step_types;
    [all_pts_f{j}, all_pts_r{j}] = get_all_pts( step_types, NITER, TransformLibrary );
    [C_eff_overlap_f{j}, C_eff_overlap_r{j}] = get_C_eff_overlap( all_pts_f{j}, all_pts_r{j} );
end
%%
clf
for j = 1:N;
    plot( 1:N,C_eff_overlap_f{j} ); hold on
    plot( 1:N,C_eff_overlap_r{j}); hold on
    set(gca,'fontweight','bold','xgrid','on'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
    title(sprintf('Forward/reverse overlap molarity for circularization of %d-mer with single inline conf',N) );
    legends = [legends, {sprintf('C_eff_overlap_f, inline at %d',j),sprintf('C_eff_overlap_r, inline at %d',j)}];
end
h = legend( legends); set(h,'interpreter','none');

%% Let's see if circularized trajectories make sense
j = N/2;
sample_circle_trajectory( step_type_sets{j}, TransformLibrary, all_pts_r{j} );

%% Scan through circularization lengths...
%N_overlap = [2:100];NITER=1000;
N_overlap = [2:20];NITER=1000;
step_types = repmat({'BB'},1,max(N_overlap));
[all_pts_f_BB, all_pts_r_BB] = get_all_pts( step_types, NITER, TransformLibrary );
[C_eff_overlap_halfway_BB,C_eff_overlap_halfway_error_BB] = get_C_eff_overlap_halfway( N_overlap, all_pts_f_BB, all_pts_r_BB );

%% Check effect of inline -- should disappear for longer lengths.
step_types_Inline = repmat({'BB'},1,max(N_overlap)); step_types_Inline{1} = 'Inline';
[all_pts_f_Inline, all_pts_r_Inline] = get_all_pts( step_types, NITER, TransformLibrary );
[C_eff_overlap_halfway_Inline,C_eff_overlap_halfway_error_Inline] = get_C_eff_overlap_halfway( N_overlap, all_pts_f_Inline, all_pts_r_Inline );

%%
clf
errorbar( N_overlap, C_eff_overlap_halfway_BB, C_eff_overlap_halfway_error_BB); hold on
errorbar( N_overlap, C_eff_overlap_halfway_Inline, C_eff_overlap_halfway_error_Inline); hold on
set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
legend( 'Any BB: overlap of forward/reverse','Single Inline link: overlap of forward/reverse ' );
set(gcf, 'PaperPositionMode','auto','color','white');



