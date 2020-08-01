function [C_eff, C_eff_err,out] = sample_motif( motif, NITER, TransformLibrary );
% out = sample_motif( motif, NITER, TransformLibrary );
%
% Wrapper around GET_C_EFF_OVERLAP
%
% Inputs
%  motif = step_types, e.g. {'BP','BB','BB','BB','BB'} to seed a triloop hairpin.
%  NITER = number of trajectories to sample for each length
%  TransformLibrary = collection of TransformSets -- one must be BB.
%
% Output
%  C_eff     = best guess for C_eff, based on average of independent estimates .
%  C_eff_err =  C_eff_err, based on stdev of independent estimates.
%  out       = struct with sub-estimates: C_eff_f, C_eff_r

NITER = 5000;
step_types = motif;
[all_pts_f, all_pts_r] = get_all_pts( step_types, NITER, TransformLibrary );
%%
tic
[C_eff_overlap_f, C_eff_overlap_r,C_eff_overlap_f_err, C_eff_overlap_r_err] = get_C_eff_overlap( all_pts_f, all_pts_r );
toc

%%
out.C_eff_f = C_eff_overlap_f;
out.C_eff_r = C_eff_overlap_r;
out.C_eff_f_err = C_eff_overlap_f_err;
out.C_eff_r_err = C_eff_overlap_r_err;

C_eff_estimates = [C_eff_overlap_f(2:end), C_eff_overlap_r(2:end)];
C_eff = mean( C_eff_estimates );
C_eff_err = std( C_eff_estimates )/sqrt(length(C_eff_estimates)-1);


%
subplot(2,1,1);
cla;
plot( C_eff_overlap_f,'linew',2 ); hold on
plot( C_eff_overlap_r,'linew',2); hold on
plot( [1 length(C_eff_overlap_f)], C_eff*[1 1], 'k','linew',0.5 );
plot( [1 length(C_eff_overlap_f)], C_eff*[1 1]+C_eff_err, 'k','linew',2 );
plot( [1 length(C_eff_overlap_f)], C_eff*[1 1]-C_eff_err, 'k','linew',2 );
set(gca,'fontweight','bold','xgrid','on'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
title(sprintf('Forward/reverse C_{eff}: %f +/- %f M', C_eff, C_eff_err) );
legend( 'forward','reverse','mean');
ylim( [0 2*max(C_eff_estimates) + 5*C_eff_err] );



%%
for n = 5:8
    subplot(2,4,n);
    sample_circle_trajectory( step_types, TransformLibrary, all_pts_r );
end
