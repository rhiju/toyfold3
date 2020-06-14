%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's try ribose rings!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load toyfold3_test.mat pdbstruct
ribose_atoms = {'C4''','C3''','C2''','C1''','O4'''};
N = length(ribose_atoms);
TransformLibrary = {};
for n = 1:N;
    n1 = mod(n  ,N)+1;
    n2 = mod(n+1,N)+1;
    n3 = mod(n+2,N)+1;
    TransformSet = get_transform_set( pdbstruct, ...
        {ribose_atoms{n}, ribose_atoms{n1},ribose_atoms{n2}},...
        {ribose_atoms{n1},ribose_atoms{n2},ribose_atoms{n3}} );
    tags{n} = sprintf('Ribose%2sprimeTo%2sprime',ribose_atoms{n}(1:2),ribose_atoms{n1}(1:2));
    TransformLibrary = setfield( TransformLibrary, tags{n}, TransformSet );
end

%%
NITER = 2000;
N = length( ribose_atoms );
for n = 1:N
    step_types = {};
    for m = 1:N; step_types{m} = tags{ mod(m+n-2,N)+1 }; end
    step_type_sets{n} = step_types;
    C_eff_ribose(n) = compute_C_eff_circular(NITER, step_types, TransformLibrary);
end
plot( C_eff_ribose); hold on
xlabel( 'Which atom ring begins at' ); ylabel( 'C_{eff}'); 
set(gca,'xtick',[1:5]);
title( ['Ribose ring closure: ',strjoin(ribose_atoms,'-') ] );

%%
n = 1;
NITER = 2000; N = length( ribose_atoms );
[all_pts_f_ribose, all_pts_r_ribose] = get_all_pts( step_type_sets{n}, NITER, TransformLibrary );
%%
clf
for m = 1:N;
    pts = all_pts_f_ribose{m}.T6; plot3( pts(:,1), pts(:,2), pts(:,3), '.' ); axis equal; axis vis3d;    hold on
end
plot3(0,0,0,'ko','markersize',10);
    
%% Show overlap
[C_eff_overlap_ribose_f, C_eff_overlap_ribose_r,C_eff_overlap_ribose_f_err, C_eff_overlap_ribose_r_err] = ...
    get_C_eff_overlap( all_pts_f_ribose, all_pts_r_ribose ); clf;
%plot( [C_eff_overlap_f; C_eff_overlap_r]' ); hold on   
errorbar( 1:N,C_eff_overlap_ribose_f,C_eff_overlap_ribose_f_err ); hold on   
errorbar( 1:N,C_eff_overlap_ribose_r,C_eff_overlap_ribose_r_err ); hold on   
plot( N,C_eff_ribose(n),'x' );  
h = legend( 'C_eff_overlap_f','C_eff_overlap_r','overlap at 0'); set(h,'interpreter','none');
set(gca,'fontweight','bold','xgrid','on'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
title(sprintf('Forward/reverse overlap molarity for circularization of %d-mer',N) );

%% sample trajectory
sample_circle_trajectory( step_type_sets{1}, TransformLibrary, all_pts_r_ribose );


