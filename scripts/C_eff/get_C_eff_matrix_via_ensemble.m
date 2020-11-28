function [C_eff,C_eff_err] = get_C_eff_matrix_via_ensemble( step_types, TransformLibrary, NITER );
% [C_eff,C_eff_err] = get_C_eff_matrix_via_ensemble( step_types, TransformLibrary, NITER );
%
% get C_eff(i,j) matrix for Watson-Crick pairing between every
% residue i and other residue j in polymer defined by step_types
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
if ~exist( 'NITER','var' ) NITER = 1000;end;

C_eff = [];
C_eff_err = [];

N = length( step_types )  + 1;
for q = 1:N; 
    all_pts{q} = struct(); 
    all_pts_BP{q} = struct(); 
end

% Note: could separate this into separate function, and then allow 
% this function to accept precomputed ensemble -- will be better
% in future when I have magic function that goes from secondary structure
% of complex RNA to ensemble.
tic
fprintf( 'Generating ensemble with %d traces...\n', NITER );
for i = 1:NITER
    y = get_random_trace(step_types, TransformLibrary, 0);
    % Could be better to use a tensor and use einstein summations later...
    %  currently doable in MATLAB but needs einsum() function.
    for q = 1:N
        all_pts{q}.t(:,i)   = y.t(:,q);
        all_pts{q}.R(:,:,i) = y.R(:,:,q);
    end
end
toc

%% Add BP helper 'appendage' to each residue in trajectory.
% These are useful later in computing C_eff for forming those BP's
tic
fprintf( 'Adding single BP to ensemble with %d traces...\n', NITER );
transforms_BP = getfield( TransformLibrary, 'BP' );
ntransforms = size(transforms_BP.t,2);
for i = 1:NITER
    for q = 1:N
        j = randi(ntransforms);        
        [all_pts_BP{q}.t(:,i), all_pts_BP{q}.R(:,:,i)] = ...
            apply_transform(all_pts{q}.t(:,i),  all_pts{q}.R(:,:,i), ...
                            transforms_BP.t(:,j),transforms_BP.R(:,:,j));
    end
end
toc

tic
for i = 1:length( step_types)
    fprintf( 'Doing row %d of %d...\n',i,length(step_types));
    for j = 1:length( step_types)
        % heuristic -- half-way point between i and j
        %  provides a legitimate coordinate system from which
        %  6D distributions for i and j can be computed and then
        %  overlap integral is super-efficient way to compute C_eff.
        %
        %              (i*)
        %               | BP
        %   1---i---q---j-------N
        %    <---------- T(q->1)
        %    ---> T(1->i)
        %
        %  T(q->i)  = T(q->1) x T(1->i) = T(1->q)[reversed] x T(1->i)
        %  T(q->i*) = T(q->1) x T(1->i*)
        %           = T(1->q)[reversed] x T(1->j) x T(j->i*)
        %           = T(1->q)[reversed] x T(1->j) x T_BP
        %
        q = floor((i+j)/2);        
        % Needed to get into q's coordinate system... T(1->q)
        [t_rev,R_rev] = reverse_transform(all_pts{q}.t,all_pts{q}.R );
        if (i < j )
            % again einsum() might read more simply and save time.
            for n = 1:NITER
                [pts_f.t(:,n), pts_f.R(:,:,n)] = apply_transform(t_rev(:,n),R_rev(:,:,n),all_pts{i}.t(:,n),all_pts{i}.R(:,:,n) );
                [pts_r.t(:,n), pts_r.R(:,:,n)] = apply_transform(t_rev(:,n),R_rev(:,:,n),all_pts_BP{j}.t(:,n),all_pts_BP{j}.R(:,:,n) );
            end
            [C_eff(i,j), C_eff_err(i,j)] = get_C_eff_from_pts_6D(pts_f,pts_r,0,1);
        else
            for n = 1:NITER
                [pts_f.t(:,n), pts_f.R(:,:,n)] = apply_transform(t_rev(:,n),R_rev(:,:,n),all_pts_BP{j}.t(:,n),all_pts_BP{j}.R(:,:,n) );
                [pts_r.t(:,n), pts_r.R(:,:,n)] = apply_transform(t_rev(:,n),R_rev(:,:,n),all_pts{i}.t(:,n),all_pts{i}.R(:,:,n) );
            end
            [C_eff(i,j), C_eff_err(i,j)] = get_C_eff_from_pts_6D(pts_f,pts_r,0,1);
        end
    end
end
toc
