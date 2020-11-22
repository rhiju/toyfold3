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

%% 
motifs = {};
% apical loops
for i= 3:12; motifs = [ motifs, {[i]}]; end
% 2WJ
for L= 0:8; 
    for i = 0:L
        j = L-i;
        motifs = [ motifs, {[i,j]}];
    end
end
% 3WJ
for L= 0:4;
    for i = 0:L
        for j = 0:L
            k = L-i-j;
            if ( k < 0 ); continue; end;
            new_motif = [i,j,k];
            % check for uniqueness
            ok = 1;
            for q = 1:length(motifs)
                if length(motifs{q})~=3; continue;end
                for n = 1:3
                    if all(circshift(new_motif,n)==motifs{q}); ok = 0; break; end;
                end
            end
            if ~ok; continue;end;
            motifs = [ motifs, {new_motif}];
        end
    end
end

% 4WJ
for L= 0:4;
    for i = 0:L
        for j = 0:L
            for k = 0:L
                l = L-i-j-k;
                if ( l < 0 ); continue; end;
                new_motif = [i,j,k,l];
                % check for uniqueness
                
                ok = 1;
                for q = 1:length(motifs)
                    if length(motifs{q})~=4; continue;end
                    for n = 1:4
                        if all(circshift(new_motif,n)==motifs{q}); ok = 0; break; end;
                    end
                end
                if ~ok; continue;end;
                motifs = [ motifs, {new_motif}];
            end
        end
    end
end
%%
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

%%
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
make_lines

