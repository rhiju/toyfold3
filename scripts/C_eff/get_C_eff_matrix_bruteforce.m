function [C_eff,C_eff_err] = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER );
% [C_eff,C_eff_err] = get_C_eff_matrix_bruteforce( step_types_all, TransformLibrary, NITER );
%
% get C_eff(i,j) matrix for Watson-Crick pairing between every
% residue i and other residue j in polymer defined by step_types. Computes
% polymer ensembles separately for each i,j. See
% GET_C_EFF_MATRIX_VIA_ENSEMBLE for more efficient calculation. (about 2x faster)
%
% Inputs
%  step_types = list of steps ('BB',etc.) (Number of nucleotides N is length of this list plus 1 )
%  TransformLibrary = collection of TransformSets with field names corresponding to
%       step_types.
%  NITER = number of ensemble members to compute [default 1000]
%
% Outputs
%  C_eff = [N x N] Effective molarity for Watson Crick pairing (in M) 
%  C_eff_err = [N X N ] error estimate in C_eff.  
%
%
% (C) R. Das, Stanford, 2020
if~exist( 'NITER','var' ) NITER = 200;end;

C_eff = [];
C_eff_err = [];

for i = 1:length( step_types_all)
    fprintf( 'Doing row %d of %d...\n',i,length(step_types_all));
    for j = 1:length( step_types_all)
        if (i < j )
            step_types = [step_types_all(i:j), {'BP'}];
            [C_eff(i,j),C_eff_err(i,j)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
        else
            step_types = [{'BP'},step_types_all(j:i)];
            [C_eff(i,j),C_eff_err(i,j)] = get_C_eff_overlap_halfway( step_types, TransformLibrary, NITER );
        end
    end
end