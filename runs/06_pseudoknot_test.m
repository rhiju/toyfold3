%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get dG for  a-b-c-d pseudoknot:
%
%       _ BP
%      /b|--|d 
%     |  |   \ 
%   a |--|c   |
%      BP|___/
%
%  a--c and b--d are base pairs
%
% * Use 1M 'standard state'; getting dG (loop entropy).
% * To get full dG need to take into account K_d's of BP's as follows:
%
%     dG = dG (loop entropy) - Sum [ kT log( 1M / K_d(BP i) ) ]
%                              BP i
%
% * get dG (in k_B T) 
% * This script will do calculation in an ad hoc way -- later can
%    generalize to allow simpler input of build-up path
% * Most important test -- different build-up paths should give
%    same answer.
% * Let's imagine a--b is 2 nt, b--c is 4 nt, c--d is 3 nt. Then
%    there should be coupling  between left hand loop formation (a-b-c) and
%      right hand loop formation (b-c-d). 
%     And.. test that buildup of left loop, then right; vs. right loop
%       then left give same dG  -- should be non-trivial!
% * Can also test paths like a-b-d-c loop first, then connect b-c.
% * Can also calculate properties of ensemble.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in library.
pdbstruct = pdbread( '../data/4ybb_DIII.pdb');
stems = read_stems_toyfold3( '../data/4ybb_DIII.pdb.stems.txt' );
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
BB_dinucleotides_reverse = reverse_dinucleotides( BB_dinucleotides );
TransformLibrary.BBr = get_transform_set( pdbstruct, BB_dinucleotides_reverse, {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% First simulate a-c-b-d as an
%  acyclic tree spanning the pseudoknot.
%         
%       b|--|d 
%        |
%        |
%        |
%   a |--|c   
%         
%  ([...)]
%
%  c-b=4 is 5 nts -- note backwards steps though!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%step_types_acbd = {'BP','BBr','BBr','BBr','BBr','BP'};
step_types_acbd = {'BP','BBr','BBr','BBr','BBr','BBr','BP'};
%NITER = 50000; clear tree_traces;
NITER = 1000; clear tree_traces;
for i = 1:NITER
    tree_traces{i} = get_random_trace(step_types_acbd, TransformLibrary, 0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let's imagine a--b is 1 nt, b--c is 5 nt, c--d is 1 nt:
%
%  ([....)] --> pseudoknot
%
% Now determine C_eff for closing each chain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
step_types_f_ab = {'BB'};
pts_f_ab = get_pts_forward( 1000, step_types_f_ab, TransformLibrary);
pts_r_ab = get_transforms_from_traces(tree_traces,1,length(step_types_acbd));
[~,~,C_eff_ab_samples] = get_C_eff_from_pts_6D( pts_r_ab, pts_f_ab );

%%
step_types_f_cd = {'BB'};
pts_f_cd = get_pts_forward( 1000, step_types_f_cd, TransformLibrary);
pts_r_cd = get_transforms_from_traces(tree_traces,2,1+length(step_types_acbd));
[~,~,C_eff_cd_samples] = get_C_eff_from_pts_6D( pts_r_cd, pts_f_cd );

%%
clf
plot( C_eff_ab_samples, C_eff_cd_samples, '.' );
xlabel( 'C_{eff}(ab)');
ylabel( 'C_{eff}(cd)');
%  (....)  --> pentaloop
fprintf( '\n\nC_eff (ab): %5.3f\n', mean( C_eff_ab_samples ))
%   [....] --> pentaloop
fprintf( 'C_eff (cd): %5.3f\n', mean( C_eff_cd_samples ))
%  ([...)] --> pseudoknot
fprintf( 'K=C_eff(ab,cd)/C_eff(ab)C_eff(cd): %5.3f\n', ...
    mean( C_eff_ab_samples.*C_eff_cd_samples )/mean(C_eff_ab_samples)/mean(C_eff_cd_samples));
%%
%[C_eff_tetraloop,C_eff_tetraloop_err] = sample_motif( {'BP','BB','BB','BB','BB','BB'}, 2000, TransformLibrary );
%%
%[C_eff_tetraloop_rev,C_eff_tetraloop_rev_err] = sample_motif( {'BP','BBr','BBr','BBr','BBr','BBr'}, 5000, TransformLibrary );

%%
[C_eff_pentaloop,C_eff_pentaloop_err] = sample_motif( {'BP','BB','BB','BB','BB','BB','BB'}, 2000, TransformLibrary );
%%
[C_eff_pentaloop_rev,C_eff_pentaloop_err] = sample_motif( {'BP','BBr','BBr','BBr','BBr','BBr','BBr'}, 5000, TransformLibrary );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hmm, maybe its the single-nt a-b and c-d spacing throwing off KDE. 
% Try longer spacings
%
%  (..[..)..] --> pseudoknot
%  a  b  c  d
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_types_acbd = {'BP','BBr','BBr','BBr','BP'};
NITER = 500000; clear tree_traces;
for i = 1:NITER
    tree_traces{i} = get_random_trace(step_types_acbd, TransformLibrary, 0);
end

step_types_f_ab = {'BB','BB','BB'};
pts_f_ab = get_pts_forward( 1000, step_types_f_ab, TransformLibrary);
pts_r_ab = get_transforms_from_traces(tree_traces,1,length(step_types_acbd));
[~,~,C_eff_ab_samples] = get_C_eff_from_pts_6D( pts_r_ab, pts_f_ab );

step_types_f_cd = {'BB','BB','BB'};
pts_f_cd = get_pts_forward( 1000, step_types_f_cd, TransformLibrary);
pts_r_cd = get_transforms_from_traces(tree_traces,2,1+length(step_types_acbd));
[~,~,C_eff_cd_samples] = get_C_eff_from_pts_6D( pts_r_cd, pts_f_cd );

clf
plot( C_eff_ab_samples, C_eff_cd_samples, '.' );
xlabel( 'C_{eff}(ab)');
ylabel( 'C_{eff}(cd)');
%  (....)  --> pentaloop
fprintf( '\n\nC_eff (ab): %5.3f\n', mean( C_eff_ab_samples ))
%   [....] --> pentaloop
fprintf( 'C_eff (cd): %5.3f\n', mean( C_eff_cd_samples ))
%  ([...)] --> pseudoknot
fprintf( 'K=C_eff(ab,cd)/C_eff(ab)C_eff(cd): %5.3f\n', ...
    mean( C_eff_ab_samples.*C_eff_cd_samples )/mean(C_eff_ab_samples)/mean(C_eff_cd_samples));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try importance sampling
%
% Again model
%  (..[..)..] --> pseudoknot
%  a  b  c  d
%
% But try instead to compare:
%
%     [.....] --> pentaloop
%     b     d
%
%  to b--c segemnt sampled from pre-circularized abc:
%
%  (..x..)
%  a  b  c
%     |__|
%       |
%       v
%      __
%     |  |
%     [.....] --> pentaloop
%     b  c  d
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_types_abc = {'BB','BB','BB','BB','BB','BB','BP'};
NITER = 500;
[all_pts_f_abc, all_pts_r_abc] = get_all_pts( step_types_abc, NITER, TransformLibrary );
%%
clear traces_circle_abc;
for i = 1:200
    traces_circle_abc{i} = sample_circle_trajectory( step_types_abc, TransformLibrary, all_pts_r_abc );
end
pts_f_bc_circle = get_transforms_from_traces(traces_circle_abc,3,6); % extract b-->c steps from within circle
%%
step_types_f_bc = {'BB','BB','BB'};
pts_f_bc = get_pts_forward( 1000, step_types_f_bc, TransformLibrary);

%%
step_types_r_bc = {'BP','BBr','BBr','BBr'};
pts_r_bc = get_pts_forward( 1000, step_types_r_bc, TransformLibrary);

%%
[C_eff_bc,C_eff_err_bc,C_eff_bc_samples] = get_C_eff_from_pts_6D( pts_f_bc, pts_r_bc );
C_eff_relerr_bc = C_eff_err_bc/C_eff_bc;

[C_eff_bc_circle,C_eff_err_bc_circle,C_eff_bc_circle_samples] = get_C_eff_from_pts_6D( pts_f_bc_circle, pts_r_bc );
C_eff_relerr_bc_circle = C_eff_err_bc_circle/C_eff_bc_circle;

C_eff_ratio = C_eff_bc_circle/C_eff_bc;
C_eff_relerr_ratio = sum( C_eff_relerr_bc^2 + C_eff_relerr_bc_circle^2 );
C_eff_err_ratio = C_eff_ratio * C_eff_relerr_ratio;

fprintf('\n');
fprintf( 'C_eff_bcd              %f+/-%f\n',C_eff_bc,C_eff_err_bc);
fprintf( 'C_eff_bcd [circle abc] %f+/-%f\n',C_eff_bc_circle,C_eff_err_bc_circle);

fprintf( 'Ratio [circle/lin]     %f+/-%f\n',C_eff_ratio,C_eff_err_ratio);
