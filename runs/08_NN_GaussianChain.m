%% Sanity checks -- check 'general' Gaussian chain function expression to 
% special cases for 1,2,3,4 rods

C_eff_gaussian_chain_func( [5 ], 1, 0 )
C_eff_gaussian_chain_func( [5 ], 1, 1 )

C_eff_gaussian_chain_func( [5 5 ], 1, 0 )
C_eff_gaussian_chain_func( [5 5 ], 1, 1 )

C_eff_gaussian_chain_func( [5 5 10 ], 1, 0 )
C_eff_gaussian_chain_func( [5 5 10 ], 1, 1 )

C_eff_gaussian_chain_func( [5 5 10 10], 1, 0 )
C_eff_gaussian_chain_func( [5 5 10 10], 1, 1 )

%% Plot some scans
x = [0.1:0.1:20]; 
subplot(3,2,1);
for i = 1:length(x); y(i) = C_eff_gaussian_chain_func( [x(i)], 1 ); end; 
plot(x,y)
xlabel( 'Rod length');
ylabel( 'C_{eff} (M)');
title( 'Single rod, Gaussian variance of 1 Å');
set(gcf, 'PaperPositionMode','auto','color','white');

subplot(3,2,2);
for i = 1:length(x); y(i) = C_eff_gaussian_chain_func( [10 x(i)], 1 ); end; 
plot(x,y)
xlabel( 'Rod length');
ylabel( 'C_{eff} (M)');
title( 'Two rods, one rod is 10 Å;\newlineGaussian variance of 1 Å');
set(gcf, 'PaperPositionMode','auto','color','white');

subplot(3,2,3);
for i = 1:length(x); y(i) = C_eff_gaussian_chain_func( [10 5 x(i)], 1 ); end; 
plot(x,y)
xlabel( 'Rod length');
ylabel( 'C_{eff} (M)');
title( 'Three rods, other rods are 10 Å, 5 Å;\newlineGaussian variance of 1 Å');
set(gcf, 'PaperPositionMode','auto','color','white');

subplot(3,2,4);
for i = 1:length(x); y(i) = C_eff_gaussian_chain_func( [10 10 10 x(i)], 1 ); end; 
plot(x,y)
xlabel( 'Rod length');
ylabel( 'C_{eff} (M)');
title( 'Four rods, other rods are 10 Å;\newlineGaussian variance of 1 Å');
set(gcf, 'PaperPositionMode','auto','color','white');

subplot(3,2,5);
for i = 1:length(x); y(i) = C_eff_gaussian_chain_func( [10 5 5 5 x(i)], 1 ); end; 
plot(x,y)
xlabel( 'Rod length');
ylabel( 'C_{eff} (M)');
title( 'Five rods, other rods are 10 Å, 5 Å, 5 Å, 5 Å;\newlineGaussian variance of 1 Å');
set(gcf, 'PaperPositionMode','auto','color','white');

subplot(3,2,6);
for i = 1:length(x); y(i) = C_eff_gaussian_chain_func( [5 5 5 5 5 x(i)], 1 ); end; 
plot(x,y)
xlabel( 'Rod length');
ylabel( 'C_{eff} (M)');
title( 'Six rods, other rods are 5 Å;\newlineGaussian variance of 1 Å');
set(gcf, 'PaperPositionMode','auto','color','white');

%% Now compute, compare to ToyFold estimates.
load 07_NN_compile motifs save_C_eff  save_motif_number save_C_eff_err motif_tags
for n = 1:length( motifs )
    motif = motifs{n};
    D = 10 * ones(1,length(motif));
    % Rosetta loop_close uses 5 Å^2 for each linkage...
    gaussian_variance = sum( motif+1 ) * 6; 
    C_eff_gcf(n) = C_eff_gaussian_chain_func( D, gaussian_variance );  
end

clear C_eff_gcf_justrods;
for n = 1:length( motifs )
    motif = motifs{n};
    D = [];
    for k = 1:length(motif )
        D = [D, [20,2*(motif(k)+1)]];
    end
    % Rosetta loop_close uses 5 Å^2 for each linkage...
    gaussian_variance = sum( motif+1 ) * 2; 
    C_eff_gcf_justrods(n) = C_eff_gaussian_chain_func( D, gaussian_variance );
end

p_init = log([10,5]);
[chisq,C_eff_gcf_fit ] = loss_gaussian_chain(p_init, motifs, save_motif_number, save_C_eff,save_C_eff_err);
chisq

p_init = log([10 5]);
p = fminsearch('loss_gaussian_chain',p_init,[],motifs, save_motif_number, save_C_eff,save_C_eff_err);
[chisq,C_eff_gcf_fit ] = loss_gaussian_chain(p, motifs, save_motif_number, save_C_eff,save_C_eff_err);
chisq

set(figure(2),'pos',[10   413   884   434]);
clf
%scatter( save_motif_number, save_C_eff,'linewidth',2 );
bubblechart( save_motif_number, save_C_eff, save_C_eff./save_C_eff_err ); hold on
plot(C_eff_gcf,'o');
plot(C_eff_gcf_fit,'o');
%plot(C_eff_gcf_justrods,'o');
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


