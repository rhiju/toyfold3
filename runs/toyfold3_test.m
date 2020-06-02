%pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
pdbstruct = pdbread( '../data/4ybb_23S.pdb');
%%
[ctr,M,chainbreak] = get_frames( pdbstruct );
%cla; draw_trace( ctr, M );

%%
[t,R] = get_transform_library(ctr, M, chainbreak);

%% Cool let's generate a random trajectory
N = 100; % number of nucleotides, steps is N-1.
[x,m] = get_random_trace(N, t, R );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's get some statistics, and estimate return probability 
% (C_eff for circular RNA)
which_N = [2:9 10:5:40, 2:9 10:5:40, 50:10:100]; NITER = 100000;
%which_N = [2:12, 2:12]; NITER = 100000;
C_eff = compute_C_eff_circular(NITER, which_N, t, R);
%%
clf; plot( which_N, C_eff,'o' );set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
title('Effective molarity for circularization' );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try overlap between forward and reverse distributions as
%  potentially quite efficient C_eff calculator.
%   ... though seeing some problems, likely because of edge effects
%     at boundaries of SO(3) rotation space (+/-pi in axis-angle).
Nmax = 10; NITER = 2000;
[all_pts_f, all_pts_r] = get_all_pts( Nmax, NITER, t, R );

%% Show overlap
N = 20;
[C_eff_overlap_f, C_eff_overlap_r] = get_C_eff_overlap( N, all_pts_f, all_pts_r ); clf;
plot( [C_eff_overlap_f; C_eff_overlap_r]' ); hold on   
plot( N,mean( C_eff(find(which_N==N)) ),'x' );  
h = legend( 'C_eff_overlap_f','C_eff_overlap_r','overlap at 0'); set(h,'interpreter','none');
set(gca,'fontweight','bold'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
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

