function [chisq,C_eff_gcf] = loss_gaussian_chain(p, motifs, save_motif_number, save_C_eff,save_C_eff_err);
% [chisq,C_eff_gcf] = loss_gaussian_chain(p, motifs, save_motif_number, save_C_eff,save_C_eff_err);
length_across_helix = exp(abs(p(1)));
gaussian_variance_per_nt = exp(abs(p(2)));
for n = 1:length( motifs )
    motif = motifs{n};
    D = length_across_helix * ones(1,length(motif));
    % Rosetta loop_close uses 5 Å^2 for each linkage...
    gaussian_variance = sum( motif+1 ) * gaussian_variance_per_nt; 
    C_eff_gcf(n) = C_eff_gaussian_chain_func( D, gaussian_variance );  
end
save_C_eff_relerr = save_C_eff_err./save_C_eff;
gp = find( ~isnan(save_C_eff) & save_C_eff_relerr>0 );

chisq = sum(  (log(C_eff_gcf(save_motif_number(gp))) - log(save_C_eff(gp)) ).^2 ./ save_C_eff_relerr(gp).^2);
%chisq = sum(  (log(C_eff_gcf(save_motif_number(gp))) - log(save_C_eff(gp)) ).^2 );