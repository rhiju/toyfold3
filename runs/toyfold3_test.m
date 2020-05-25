%pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
pdbstruct = pdbread( '../data/4ybb_23S.pdb');
%%
[ctr,M,chainbreak] = get_frames( pdbstruct );
cla; draw_trace( ctr, M );

%%
[t,R] = get_transform_library(ctr, M, chainbreak);

%% Cool let's generate a random trajectory
N = 100;
[x,m] = get_random_trace(N, t, R );

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
    

    
    
    
