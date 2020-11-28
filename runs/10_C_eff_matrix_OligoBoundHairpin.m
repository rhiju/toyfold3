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

%%  
% New ... needed to construct helices.
BB_stem_dinucleotides = get_BB_from_stems( stems );
TransformLibrary.BB_stem = get_transform_set( pdbstruct, BB_stem_dinucleotides,  {'C5''','C4''','C3'''},{'C5''','C4''','C3'''} );
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualize transforms -- should be pretty regular
clf
plot( TransformLibrary.BB_stem.T6(:,2),TransformLibrary.BB_stem.T6(:,5),'o'); hold on
plot( TransformLibrary.BB.T6(:,2),TransformLibrary.BB.T6(:,5),'.'); 
xlabel( 'y'); ylabel( 'v_y');
legend( 'BB_{stem}', 'BB');

%% Quick test of long helix traces.
step_types = repmat( {'BB_stem'},1,20); % number of nucleotides, steps is N-1.
x = get_random_trace(step_types, TransformLibrary ); 

%% C_eff matrix for
%
%        GGGGG
%        |||||
% ...AAAACCCCCAAAUU...
% |___________|
%    C_eff(i,j)
NITER = 500;
step_types_all = [repmat({'BB'},1,25)];
C_eff_25mer = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER);
%%
step_types_all = [repmat({'BB'},1,10),repmat({'BB_stem'},1,5),repmat({'BB'},1,10)];
C_eff_25mer_oligo5bound = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER);

%%
step_types_all = [repmat({'BB'},1,30)];
C_eff_30mer = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER);
step_types_all = [repmat({'BB'},1,10),repmat({'BB_stem'},1,10),repmat({'BB'},1,10)];
C_eff_30mer_oligo10bound = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER);

%%
set(figure(1),'pos',[13   454   367   350]);
subplot(2,2,1);
imagesc( C_eff_25mer,[0 0.05] )
make_lines_horizontal([0 10 15 25]);
make_lines([0 10 15 25]);
title( '25mer')

subplot(2,2,2);
imagesc( C_eff_25mer_oligo5bound,[0 0.05] )
make_lines_horizontal([0 10 15 25]); make_lines([0 10 15 25]);
title( '25mer oligo5 bound')

subplot(2,2,3);
imagesc( C_eff_30mer,[0 0.05] )
make_lines_horizontal([0:10:30]);
make_lines([0:10:30]);
title( '30mer')

subplot(2,2,4);
imagesc( C_eff_30mer_oligo10bound,[0 0.05] )
make_lines_horizontal([0:10:30]);
make_lines([0:10:30]);
title( '30mer oligo10 bound')
set(gcf, 'PaperPositionMode','auto','color','white');

%% Infer from ensemble

%% Pilot
NITER = 500;
step_types_all = [repmat({'BB'},1,10)];
tic
C_eff_10mer = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER);
toc
%%
tic
C_eff_10mer_TEST = get_C_eff_matrix_via_ensemble( step_types_all, TransformLibrary, NITER);
toc

%%
subplot(1,2,1);
imagesc( C_eff_10mer,[0 0.05] );
title( '10mer -- bruteforce');
subplot(1,2,2);
imagesc( C_eff_10mer_TEST,[0 0.05] );
title( '10mer -- reuse ensemble');
set(gcf, 'PaperPositionMode','auto','color','white');



%%
NITER = 2000;
step_types_all = [repmat({'BB'},1,10),repmat({'BB_stem'},1,5),repmat({'BB'},1,10)];
C_eff_25mer_oligo5bound_TEST = get_C_eff_matrix_via_ensemble( step_types_all, TransformLibrary, NITER);

%%
NITER = 2000;
step_types_all = [repmat({'BB'},1,10),repmat({'BB_stem'},1,5),repmat({'BB'},1,10)];
C_eff_25mer_oligo5bound_bruteforce2k = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER);

%%
set(figure(2),'pos',[57   326   570   228]);
subplot(1,2,1);
imagesc( C_eff_25mer_oligo5bound_bruteforce2k,[0 0.05] );
title( '25mer oligo5bound -- bruteforce');
make_lines_horizontal([0 10 15 25]); make_lines([0 10 15 25]);

subplot(1,2,2);
imagesc( C_eff_25mer_oligo5bound_TEST,[0 0.05] );
title( '25mer oligo5bound  -- reuse ensemble');
make_lines_horizontal([0 10 15 25]); make_lines([0 10 15 25]);
set(gcf, 'PaperPositionMode','auto','color','white');


%%
C_eff_25mer_oligo5bound_ensemble2k = C_eff_25mer_oligo5bound_TEST;
NITER = 2000;
step_types_all = repmat({'BB'},1,25);
C_eff_25mer_ensemble2k = get_C_eff_matrix_via_ensemble( step_types_all, TransformLibrary, NITER);

% OK now, let's compute log-ratio with and without oligo bound!
%% log-odds --> legit calculation for 25mer
set(figure(3),'pos',[57   126   570   228]);
subplot(1,2,1);
imagesc( log(C_eff_25mer_ensemble2k),[-6 -3] );
title( '25mer -- ensemble2k');
make_lines_horizontal([0 10 15 25]); make_lines([0 10 15 25]);

subplot(1,2,2);
imagesc( log(C_eff_25mer_oligo5bound_ensemble2k),[-6 -3] );
title( '25mer oligo5bound -- ensemble2k');
make_lines_horizontal([0 10 15 25]); make_lines([0 10 15 25]);
set(gcf, 'PaperPositionMode','auto','color','white');

%%
set(figure(4),'pos',[ 0   109   277   225]);
imagesc( -log(C_eff_25mer_oligo5bound_ensemble2k./C_eff_25mer_ensemble2k),[-2 2] );
title( '\Delta\DeltaG = -log(oligo-bound/no-oligo))');
set(gcf, 'PaperPositionMode','auto','color','white');
make_lines_horizontal([0 10 15 25]); make_lines([0 10 15 25]);
colormap( customcolormap_preset('red-white-blue') );
colorbar();

%% 
% OK now, let's compute log-ratio with and without oligo bound!
NITER = 2000;
step_types_50mer = repmat({'BB'},1,50);
C_eff_50mer_ensemble2k = get_C_eff_matrix_via_ensemble( step_types_50mer, TransformLibrary, NITER);

%%
step_types_50mer_oligo10bound = step_types_50mer;
for i = 21:30; step_types_50mer_oligo10bound{i} = 'BB_stem';end;
C_eff_50mer_oligo10bound_ensemble2k = get_C_eff_matrix_via_ensemble( step_types_50mer_oligo10bound, TransformLibrary, NITER);

%%
set(figure(3),'pos',[57   126   570   228]);
subplot(1,2,1);
imagesc( log(C_eff_50mer_ensemble2k),[-6 -3] );
title( 'C_{eff} 50mer');
bounds = [0 20 30 50];
make_lines_horizontal(bounds); make_lines(bounds);

subplot(1,2,2);
imagesc( log(C_eff_50mer_oligo10bound_ensemble2k),[-6 -3] );
title( 'C_{eff} 50mer [10-mer double helix in center]');
make_lines_horizontal(bounds); make_lines(bounds);
set(gcf, 'PaperPositionMode','auto','color','white');

%%
set(figure(4),'pos',[ 0   139   327   280]);
imagesc( -log(C_eff_50mer_oligo10bound_ensemble2k./C_eff_50mer_ensemble2k),[-2 2] );
title( '\Delta\DeltaG [50mer, oligo10] = -log(oligo-bound/no-oligo))');
set(gcf, 'PaperPositionMode','auto','color','white');
make_lines_horizontal(bounds); make_lines(bounds);
colormap( customcolormap_preset('red-white-blue') );
colorbar();

