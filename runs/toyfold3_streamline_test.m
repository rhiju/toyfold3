%%
pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
%%
tic
TransformSet = get_transform_set( pdbstruct, {'C5''','C4''','C3'''},{'C5''','C4''','C3''',1} );
toc

%%
TransformLibrary = {};
TransformLibrary = setfield( TransformLibrary, 'BB', TransformSet );

%% For ribose rings
%TransformSet = get_transform_library( pdbstruct, {'O4''','C4''','C3'''},{'C4''','C3''','C2'''} );
%Tranformlibrary = add_transform_library( Tranformlibrary, TransformSet, 'RiboseC4primeToC3prime' ); % backbone

%% for stems
% base_pairs = get_base_pairs('4ybb_DIII.stems.txt');
% TransformSet = get_transform_library( pdbstruct,  {'C5''','C4''','C3'''},{'C5''','C4''','C3''','base_pair'}, base_pairs );

%% Cool let's generate a random trajectory
step_types = repmat( {'BB'},1,100); % number of nucleotides, steps is N-1.
x = get_random_trace(step_types, TransformLibrary );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's get some statistics, and estimate return probability 
% (C_eff for circular RNA)
%which_N = [2:19 20:5:40, 2:19 20:5:40, 2:19, 2:19, 50:10:100]; NITER = 5000;
%which_N = [2:12, 2:12]; NITER = 100000;
%C_eff = compute_C_eff_circular(NITER, which_N, t, R);
%%
%clf; plot( which_N, C_eff,'o' );set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
%title('Effective molarity for circularization' );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try overlap between forward and reverse distributions as
%  potentially quite efficient C_eff calculator.
%   ... though seeing some problems, likely because of edge effects
%     at boundaries of SO(3) rotation space (+/-pi in axis-angle).
Nmax = 10; NITER = 1000;
[all_pts_f, all_pts_r] = get_all_pts( Nmax, NITER, TransformLibrary );

%% Show overlap
N = 10;
[C_eff_overlap_f, C_eff_overlap_r] = get_C_eff_overlap( N, all_pts_f, all_pts_r ); clf;
plot( [C_eff_overlap_f; C_eff_overlap_r]' ); hold on   
plot( N,mean( C_eff(find(which_N==N)) ),'x' );  
h = legend( 'C_eff_overlap_f','C_eff_overlap_r','overlap at 0'); set(h,'interpreter','none');
set(gca,'fontweight','bold','xgrid','on'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
title(sprintf('Forward/reverse overlap molarity for circularization of %d-mer',N) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try pure overlap based calc as function of Nsteps
N_overlap = [2:100];
C_eff_overlap_halfway = get_C_eff_overlap_halfway( N_overlap, all_pts_f, all_pts_r );

%%
clf
plot( which_N, C_eff,'o' ); hold on
plot( N_overlap, C_eff_overlap_halfway); hold on
set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
legend( 'circularize to origin','overlap of forward/reverse' );
set(gcf, 'PaperPositionMode','auto','color','white');

