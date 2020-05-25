%pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
pdbstruct = pdbread( '../data/4ybb_23S.pdb');
%%
[ctr,M,chainbreak] = get_frames( pdbstruct );
%cla; draw_trace( ctr, M );

%%
[t,R] = get_transform_library(ctr, M, chainbreak);

%% Cool let's generate a random trajectory
N = 100;
[x,m] = get_random_trace(N, t, R );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's get some statistics, and estimate return probability 
% (C_eff for circular RNA)
which_N = [2:9 10:5:40, 2:9 10:5:40, 50:10:100];
NITER = 20000;
clear p
for i = 1:length(which_N)
    N = which_N(i); % number of steps
    tic
    pts_f = get_pts_forward( NITER, N, t, R);
    toc
    %cla;plot3( pts_f(1,:),pts_f(2,:),pts_f(3,:),'o'); axis equal
    % then collect histograms
        
    s = get_kde_bandwidth( pts_f );
    pts_r = [0,0,0,0,0,0];

    p(i) =  mvksdensity(pts_f,pts_r,'Bandwidth',s)';
end

%%
clf
C_eff = p/(1/(8*pi^2)*6.022e23/1e27 )
plot( which_N, C_eff,'o' );set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
title('Effective molarity for circularization' );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try overlap between forward and reverse distributions as
%%  potentially quite efficient C_eff calculator.
Nmax = 100;
NITER = 2000;
for i = 88:Nmax
    tic
    all_pts_f{i} = get_pts_forward( NITER, i, t, R);
    % don't really need to repeat above -- could just reverse transforms
    all_pts_r{i} = get_pts_reverse( NITER, i, t, R); 
    toc
end

%%
clear C_eff_overlap_f C_eff_overlap_r
N = 10;
for i = 1:N-1
    C_eff_overlap_f(i) = get_C_eff_from_pts_6D(all_pts_f{i},all_pts_r{N-i});
    C_eff_overlap_r(i) = get_C_eff_from_pts_6D(all_pts_r{i},all_pts_f{N-i});
end
C_eff_overlap_f(N) = get_C_eff_from_pts_6D(all_pts_f{i},[0,0,0,0,0,0]);
C_eff_overlap_r(N) = get_C_eff_from_pts_6D(all_pts_r{i},[0,0,0,0,0,0]);
%%
clf;
N = 10;
plot( C_eff_overlap_f ); hold on   
plot( C_eff_overlap_r )    
plot( N,mean( C_eff(find(which_N==N)) ),'x' );  
h = legend( 'C_eff_overlap_f','C_eff_overlap_r','overlap at 0'); set(h,'interpreter','none');
set(gca,'fontweight','bold'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
title(sprintf('Forward/reverse overlap molarity for circularization of %d-mer',N) );


%% try pure overlap based calc as function of Nsteps
N_overlap = [2:100];
clear C_eff_overlap_mean;
for k = 1:length(N_overlap)
    N = N_overlap(k);
    i = floor(N/2);
    C_eff_overlap_mean(k) = get_C_eff_from_pts_6D(all_pts_f{i},all_pts_r{N-i});
end
%%
clf
plot( which_N, C_eff,'o' ); hold on
plot( N_overlap, C_eff_overlap_mean); hold on
set(gca,'fontweight','bold');xlabel('N');ylabel('C_{eff} (M)');
legend( 'circularize to origin','overlap of forward/reverse' );

