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
%% simple tests first..
%% Let's take a look at closed trajectories
[C_eff_tetraloop,C_eff_tetraloop] = sample_motif( {'BP','BB','BB','BB','BB','BB'}, 1000, TransformLibrary );

%% Wrapper around C_eff_overlap_halfway
loop_lengths = [1:10];
NITER = 2000;
out_HP = scan_loop_length(loop_lengths,NITER,TransformLibrary,{'BP','BB'});
clf
errorbar( loop_lengths, out_HP.C_eff,  out_HP.C_eff_err,'linew',2); hold on

%% Let's do it! Enumerate motifs to compute C_eff for.
motifs = {};
% apical loops
for i= 1:20; motifs = [ motifs, {[i]}]; end
% 2WJ
for L= 0:8; 
    for i = 0:L
        j = L-i;
        motifs = [ motifs, {[i,j]}];
    end
end
motifs = get_unique_motifs(3, 4, motifs); % 3WJ
motifs = get_unique_motifs(4, 4, motifs); % 4WJ
motifs = get_unique_motifs(5, 4, motifs); % 5WJ
motifs = get_unique_motifs(6, 4, motifs); % 6WJ

%% Do the computation
save_C_eff       = [];
save_C_eff_err   = [];
save_motif_number = [];
save_motif_tag  = {};
save_step_types = {};
motif_tags = {};
NITER = 2000;
fprintf('\n\n');
for n = 1:length(motifs)
    motif =  motifs{n};
    motif_tag = strrep(int2str(motif),'  ','-');
    motif_tags = [motif_tags, {motif_tag}];
    step_types = {'BP'};
    for i = 1:length(motif);
        step_types = [step_types,repmat({'BB'},1,1+motif(i))];
        if (i < length(motif)) step_types = [step_types, {'BP'}]; end
    end
    fprintf( 1, 'Doing motif %s, %d of %d...\n',motif_tag,n,length(motifs));
    for q = 0:length(step_types)-1;
        step_types_circshift = circshift( step_types, q );
        [C_eff,C_eff_err] = get_C_eff_overlap_halfway( step_types_circshift, TransformLibrary, NITER );

        save_C_eff = [save_C_eff, C_eff];
        save_C_eff_err = [save_C_eff_err, C_eff_err];
        save_motif_number = [save_motif_number,n];
        save_motif_tag = [save_motif_tag,{motif_tag}];
        save_step_types = [save_step_types,{step_types_circshift}];
    end
end

%% Plot as bubble chart.
set(figure(2),'pos',[10   413   884   434]);
clf
%scatter( save_motif_number, save_C_eff,'linewidth',2 );
bubblechart( save_motif_number, save_C_eff, save_C_eff./save_C_eff_err );
xlabel( 'Motif')
ylabel( 'C_{eff} (M)');
bubblelim([1 1000])
title( 'Nearest neighbor circularization C_{eff}');
xlim([0+0.5 length(motifs)+0.5]);
set(gca,'xtick', [1:length(motif_tags)],'xticklabel',motif_tags,'xticklabelrot',90 )
set(gca,'yscale','log','fontweight','bold');
set(gcf, 'PaperPositionMode','auto','color','white');

save_motif_lengths = cellfun(@length,motifs);
make_lines( find( save_motif_lengths(1:end-1)~= save_motif_lengths(2:end)),'k',2 );
make_lines;

%% Output to .csv file
filename = 'NN_C_eff_toyfold3_ALL.csv';
fid = fopen( filename, 'w');
fprintf( fid, '%s,%s,%s\n','motif','C_eff','C_eff_err');
for i = 1:length(save_C_eff)
    fprintf(fid,'%s,%f,%f\n',motif_tags{save_motif_number(i)},save_C_eff(i),save_C_eff_err(i));
end
fclose(fid);
fprintf( 'Created: %s\n', filename );

%%
unique_motif_tags = unique(motif_tags,'stable');
save_C_eff_relerr = save_C_eff_err./save_C_eff;
for i = 1:length( unique_motif_tags )
    idx = find(strcmp( motif_tags(save_motif_number), unique_motif_tags{i} ) );
    log_C_eff_mean = sum( log(save_C_eff(idx))./save_C_eff_relerr(idx).^2) ./ sum(1./save_C_eff_relerr(idx).^2);
    C_eff_mean_relerr = sqrt( 1/ sum(1./save_C_eff_relerr(idx).^2) );
    C_eff_mean(i) = exp( log_C_eff_mean );
    C_eff_mean_err(i) = C_eff_mean(i) * C_eff_mean_relerr;
end

filename = 'NN_C_eff_toyfold3.csv';
fid = fopen( filename, 'w');
fprintf( fid, '%s,%s,%s\n','motif','C_eff','C_eff_err');
for i = 1:length(unique_motif_tags)
    fprintf(fid,'%s,%f,%f\n',unique_motif_tags{i},C_eff_mean(i),C_eff_mean_err(i));
end
fclose(fid);
fprintf( 'Created: %s\n', filename );





