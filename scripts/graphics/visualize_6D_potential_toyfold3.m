function visualize_6D_potential( T, boxsize, apply_neg_log )
% Shows slices of 6D potential as 2D (vx,vy) slices splayed out with
%   centers at (x,y). 
%  (Assumes z = 0, vz = 0, though perhaps we could
%    iterate over at least z stacks.)
%
% INPUTS
%   T       = struct containing tensor (6D tensor) and json with histogram info
%   boxsize = boxsize to use for visualization (default: auto)
%   apply_neg_log = show -log( p )  (default: True )
%

if ~exist( 'apply_neg_log', 'var' ) apply_neg_log = 1; end;

h_size = size( T.tensor );
xbinwidth = T.json.binwidth(1);
xbins = [T.json.minval(1) : xbinwidth : T.json.maxval(1) ];

if ~exist( 'boxsize','var' ) | boxsize == 0 | isempty( boxsize)
    % figure out longest distance from origin at which there are counts. 
    [X,Y,Z] = ndgrid( xbins, xbins, xbins );
    R = sqrt( X.^2 + Y.^2 + Z.^2 );
    if apply_neg_log; P = T.tensor; else; P = exp( -T.tensor ); end;
    P( isnan( P ) ) = 0.0;
    P = P - min(P(:)); % in case we applied a ceiling.
    hxyz = sum(sum(sum( P, 4),5),6); % sum over rotations
    boxsize = max( R( find( hxyz > 0 ) ) );
    boxsize = xbinwidth*ceil(boxsize/xbinwidth);
    boxsize = min(boxsize,max(xbins));
    fprintf( 'Visualizing with boxsize of %f, compared to histogram boxsize of %f\n', boxsize, max(xbins) );
end;

b0 = interp1( xbins, 1:h_size(1), 0, 'nearest' );
b1 = interp1( xbins, 1:h_size(1), -boxsize, 'nearest' );
b2 = interp1( xbins, 1:h_size(1), +boxsize, 'nearest' );
bins = [b1:b2];

vxbins = [T.json.minval(4) : T.json.binwidth(4) : T.json.maxval(4) ];
vb0 = interp1( vxbins, 1:h_size(4), 0, 'nearest' );
Nv = length( vxbins );

% this is pretty manual -- probably a super-slick way to do it with
% tensor manipulations
imagex = [];

Tplot = T.tensor;
if apply_neg_log; Tplot = -log( T.tensor ); end;
for i = 1:length(bins)
    xoff = Nv*(i-1);
    for j = 1:length(bins)
        yoff = Nv*(j-1);
        imagex(xoff+[1:Nv],yoff+[1:Nv]) = squeeze( Tplot(bins(i),bins(j),b0,:,:,vb0) );
    end
end
xg = interp1( 1:h_size(1), xbins, [(b1-0.5):(b2+0.5)], 'linear', 'extrap' );
xax = [min(xg)+0.5*xbinwidth/Nv ...
       max(xg)-0.5*xbinwidth/Nv];

clim = [ min( Tplot( ~isinf(Tplot) ) ) max( Tplot( ~isinf(Tplot) ) )];
fprintf( 'Min & max potential: %f %f\n', clim(1), clim(2) );
imagesc( xax, xax, imagex', clim );
hold on;
plot( 0, 0, 'kx' ); % mark origin
rotbinmax = T.json.maxval(4);
for i = 1:length(bins)
    for j = 1:length(bins)
     plot( [min(xg) max(xg)], xg(i)*[1 1],'k' );
     plot( xg(i)*[1 1], [min(xg) max(xg)],'k' );
     plot( xg(i)+xbinwidth/2, xg(j)+xbinwidth/2, 'k.' );
     xcirc = xg(i)+xbinwidth/2-(xbinwidth/2)*(180.0/rotbinmax);
     ycirc = xg(j)+xbinwidth/2-(xbinwidth/2)*(180.0/rotbinmax);
     rectangle('Position',[xcirc,ycirc,xbinwidth*(180.0/rotbinmax),xbinwidth*(180.0/rotbinmax)],'Curvature',[1 1],'EdgeColor','k');
    end
end
% let's also show vx,vy axes
i = 1; j = 1; 
xp0 = xg(i)+xbinwidth/2;
yp0 = xg(j)+xbinwidth/2;
xp = xg(i)+xbinwidth/2;
yp = xg(j)+xbinwidth/2-(xbinwidth/2)*(180.0/rotbinmax);
text( xp, yp,'\pi','VerticalAlignment','middle','HorizontalAlignment','left','Color','w');
plot( [xp0 xp],[yp0 yp], 'w-' );
xp = xg(i)+xbinwidth/2+(xbinwidth/2)*(180.0/rotbinmax);
yp = xg(j)+xbinwidth/2;
text( xp, yp,'\pi','HorizontalAlignment','center','VerticalAlignment','bottom','Color','w');
plot( [xp0 xp],[yp0 yp], 'w-' );
xp = xg(i)+xbinwidth/2;
yp = xg(j)+xbinwidth/2-(1/2)*(xbinwidth/2)*(180.0/rotbinmax);
text( xp, yp,'v_y','VerticalAlignment','middle','HorizontalAlignment','right','Color','w');
xp = xg(i)+xbinwidth/2+(1/2)*(xbinwidth/2)*(180.0/rotbinmax);
yp = xg(j)+xbinwidth/2;
text( xp, yp,'v_x','HorizontalAlignment','center','VerticalAlignment','top','Color','w');

hold off;
xt = interp1( 1:h_size(1), xbins, bins);
set(gca,'xtick',xt,'ytick',xt);
xlabel( 'x' );
ylabel( 'y' );
c = colorbar;
if apply_neg_log
    c.Label.String = '-log(p)';
else
    c.Label.String = 'E';
end
colormap( jet )
set(gcf, 'PaperPositionMode','auto','color','white');

