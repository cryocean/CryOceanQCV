function plot_map(fn,xi,yi,zi,cmax,cmin,dc,unitLab,titl, ...
    pPos,pos,interpol,lonLim,latLim)

% load coastline data
path4 = '/noc/users/cryo/QCV_Cryo2/code/';
load([path4 'gshhs_c.mat'])
xh = gshhs_c.lon;
yh = gshhs_c.lat;
lev = gshhs_c.level;
kl = (lev == 1 | lev == 6);
xh = xh(kl);
yh = yh(kl);
xh = cell2mat(xh);
yh = cell2mat(yh);
xh = xh(:)';
yh = yh(:)';

labx = sum(abs(lonLim));
laby = sum(abs(latLim));
dx = fix(labx)/6;
dy = fix(laby)/6;
DX = [2 5 10 20 30 60];
dx = abs(DX-dx);
[~,dx] = min(dx);
dy = abs(DX-dy);
[~,dy] = min(dy);
dx = DX(dx);
dy = DX(dy);
Mlab = -360:dx:360;
Plab = -90:dy:90;
    
figure(fn)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo','mapLonLimit',lonLim,'mapLatLimit',latLim);
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

% ----------------------------- with scatter ------------------------------
colormap(linspecer);
if interpol == 0
    scatterm(ha,yi,xi,6,zi,'o');
else
    hpc = pcolorm(yi,xi,zi); %#ok
end
% -------------------------------------------------------------------------
framem
mlabel('MLabelParallel','south','MLabelLocation',Mlab)
plabel('PLabelLocation',Plab)
setm(ha,'FontSize',10,'MLineLocation',60,'PLineLocation',30,'GLineWidth',4, ...
    'GAltitude',1e+30,'GLineStyle','-')
tightmap
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([cmin cmax])
set(hc,'YTick',cmin:dc:cmax)
set(hc,'FontSize',10)
ylabel(hc,unitLab) % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);

% -nodisplay mode shifts the colorbar (it is a bug) when printing. Here we
% shift back the colorbar by an amount similar to that done by Matlab
cpos(2) = cpos(2)+0.765; 
set(hc,'Position',cpos)
set(gca,'position',xpos)
gridm off
set(ha,'layer','top','TickDir','out')
title(titl,'FontSize',10)
