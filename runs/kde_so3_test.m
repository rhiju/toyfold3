%%
N = 1000;
M = get_uniform_random_matrices( N );
pts_EV = SpinCalc( 'DCMtoEV', M, 1.0e-6, 0 );
pts = [];
% convert to axis vector v_x,v_y,v_z; with length equal to rotation angle in radians
for n = 1:size( pts_EV,1); pts(n,:) = pts_EV(n,1:3) * pts_EV(n,4) * pi/180.0; end;

subplot(2,1,1);
pts_test = zeros(100,3);
vbins = [0:0.01:0.99] * pi;;
pts_test(:,3) = vbins;

s = get_kde_bandwidth( pts );
p = mvksdensity( pts, pts_test,'Bandwidth',s); 
cla
set(gcf, 'PaperPositionMode','auto','color','white');
plot( vbins, p, 'rx'); hold on
plot( vbins, (sin(vbins/2).^2)./(vbins/2).^2 * 1/(8*pi^2) );
legend( 'numerical (KDE)', '1/(8\pi^2)sinc^2(v/2)');
xlabel( 'v (angle of rotation in radians)' );
ylabel( 'Density in rotation vector space' );
ylim([0 2/(8*pi^2)]);
title( sprintf('%d samples',N) );

subplot(2,1,2);

%Apply b transform
pts_b = [];
% convert to axis vector v_x,v_y,v_z; with length equal to rotation angle in radians
v = pts_EV(:,4) * pi/180.0;
b = ( 6 * (v-sin(v))).^(1/3);
for n = 1:size( pts_EV,1); pts_b(n,:) = pts_EV(n,1:3) * b(n); end;

pts_b_test = zeros(100,3);
bbins = [0:0.01:0.99] * (6*pi)^(1/3);;
pts_b_test(:,3) = vbins;

s = get_kde_bandwidth( pts_b );
p = mvksdensity( pts_b, pts_b_test,'Bandwidth',s); 

cla
set(gcf, 'PaperPositionMode','auto','color','white');
plot( bbins, p, 'rx'); hold on
plot( bbins, 0*vbins+1/(8*pi^2) );
legend( 'numerical (KDE)', '1/(8\pi^2)');
xlabel( 'v (angle of rotation in radians)' );
ylabel( 'Density in b-vector space' );
ylim([0 2/(8*pi^2)]);

