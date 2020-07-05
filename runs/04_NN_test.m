%%
pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
%%
tic
TransformLibary = struct();
BB_dinucleotides = get_BB_dinucleotides(pdbstruct);
TransformLibrary.BB = get_transform_set( pdbstruct, BB_dinucleotides, {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc
%%
tic
stems = read_stems_toyfold3( '../data/4ybb_DIII.pdb.stems.txt' );
base_pairs = get_base_pairs_from_stems_toyfold3( stems );
TransformLibrary.BP = get_transform_set( pdbstruct, base_pairs,  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc

%% Cool let's generate a random trajectory
%% base pair step, a.k.a. stacked pair  [no closure]
step_types = {'BB','BP','BB','BP'}; 
x = get_random_trace(step_types, TransformLibrary );

%% 1x1 internal loop [no closure]
step_types = {'BB','BB','BP','BB','BB','BP'}; 
x = get_random_trace(step_types, TransformLibrary );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's get some statistics, and estimate return probability 
% (C_eff for circular RNA)
which_N = [2:19 20:5:40, 2:19 20:5:40, 2:19, 2:19, 50:10:100]; NITER = 2000;
which_N = [2:12]; NITER = 1000;
C_eff = compute_C_eff_circular_backbone(NITER, which_N, TransformLibrary);
%%
clf; plot( which_N, C_eff,'o' );set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
title('Effective molarity for circularization' );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try overlap between forward and reverse distributions as
%  potentially quite efficient C_eff calculator.
%   ... though seeing some problems, likely because of edge effects
%     at boundaries of SO(3) rotation space (+/-pi in axis-angle).
Nmax = 10; NITER = 1000;
[all_pts_f, all_pts_r] = get_all_pts( Nmax, NITER, TransformLibrary );

%% Show overlap
N = 10;
[C_eff_overlap_f, C_eff_overlap_r,C_eff_overlap_f_err, C_eff_overlap_r_err] = ...
    get_C_eff_overlap( all_pts_f(1:N), all_pts_r(1:N) ); clf;
%plot( [C_eff_overlap_f; C_eff_overlap_r]' ); hold on   
errorbar( 1:N,C_eff_overlap_f,C_eff_overlap_f_err ); hold on   
errorbar( 1:N,C_eff_overlap_r,C_eff_overlap_r_err ); hold on   
plot( N,mean( C_eff(find(which_N==N)) ),'x' );  
h = legend( 'C_eff_overlap_f','C_eff_overlap_r','overlap at 0'); set(h,'interpreter','none');
set(gca,'fontweight','bold','xgrid','on'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
title(sprintf('Forward/reverse overlap molarity for circularization of %d-mer',N) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try pure overlap based calc as function of Nsteps
N_overlap = [2:100];
[C_eff_overlap_halfway,C_eff_overlap_halfway_error] = get_C_eff_overlap_halfway( N_overlap, all_pts_f, all_pts_r );

%%
clf
plot( which_N, C_eff,'o' ); hold on
%plot( N_overlap, C_eff_overlap_halfway); hold on
errorbar( N_overlap, C_eff_overlap_halfway, C_eff_overlap_halfway_error); hold on
set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
legend( 'circularize to origin','overlap of forward/reverse' );
set(gcf, 'PaperPositionMode','auto','color','white');


%% sample trajectory
N = 6;
step_types = repmat( {'BB'},1,N); % number of nucleotides, steps is N-1. 
sample_circle_trajectory( step_types, TransformLibrary, all_pts_r );


