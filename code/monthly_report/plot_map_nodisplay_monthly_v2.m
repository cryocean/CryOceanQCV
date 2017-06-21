function plot_map_nodisplay_monthly_v2(fn,xi,yi,zi,cmax,cmin,dc,unitLab,titl, ...
    pPos,pos,xover)

if nargin < 12
    xover = false;
end
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

figure(fn)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

% ----------------------------- with scatter ------------------------------
colormap(linspecer);
if xover
%     scatterm(ha,yi,xi,10,zi,'s','filled');
    scatterm(ha,yi,xi,5,zi,'s','filled');% CHANGED 7 APRIL

else
    
    % changed for TDS
    fastscatterm_NOC(yi,xi,zi)
    
%     scatterm(ha,yi,xi,1,zi,'o');
end
% -------------------------------------------------------------------------
framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',10)
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
% % % % % % % % % % % % % % % % % % % % CHANGED FOR TDS
% % % % % % % % % % % % % % % % % % % cpos(3) = 0.6.*cpos(3);
% % % % % % % % % % % % % % % % % % % % -nodisplay mode shifts the colorbar (it is a bug) when printing. Here we
% % % % % % % % % % % % % % % % % % % % shift back the colorbar by an amount similar to that done by Matlab
% % % % % % % % % % % % % % % % % % % cpos(2) = cpos(2)+0.765; 
% % % % % % % % % % % % % % % % % % % set(hc,'Position',cpos)
set(gca,'position',xpos,'layer','top')
set(ha,'layer','top')
title(titl,'FontSize',10)
