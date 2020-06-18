% Try to do projection for a 'pseudochain' that lives only in rotation
% space SO(3)

%% The rotation part of the SE(3) rotation/translation transform 
% is just m' = R*m, same SO(3). So I can use prebuilt chains or transforms
% for SE(3) if I just pull out the rotation vectors. 
%load toyfold3_test.mat t R all_pts_f all_pts_r

Nmax = 20; NITER = 1000;
[all_pts_f, all_pts_r] = get_all_pts( Nmax, NITER, t, R );

%% Show overlap
figure(2)
N = 20;
just_SO3 = 1;
[C_eff_overlap_f, C_eff_overlap_r] = get_C_eff_overlap( N, all_pts_f, all_pts_r, just_SO3 ); clf;
plot( [C_eff_overlap_f; C_eff_overlap_r]' ); hold on   
%plot( N,mean( C_eff(find(which_N==N)) ),'x' );  
h = legend( 'C_eff_overlap_f','C_eff_overlap_r','overlap at 0'); set(h,'interpreter','none');
set(gca,'fontweight','bold'); xlabel( 'n steps for overlap'); ylabel('C_{eff} (M)');
title(sprintf('Forward/reverse overlap molarity for circularization of %d-mer',N) );
set(gcf, 'PaperPositionMode','auto','color','white');

%% take a look at the points.
figure(3)
pts = all_pts_f{19}(:,4:6);
plot3( pts(:,1), pts(:,2), pts(:,3), '.' ); axis equal; axis vis3d
hold on
pts = all_pts_r{1}(:,4:6);
plot3( pts(:,1), pts(:,2), pts(:,3), '.' ); axis equal; axis vis3d
legend( 'F 19-mer','R 1-mer');