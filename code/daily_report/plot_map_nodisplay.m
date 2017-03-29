function plot_map_nodisplay(fn,xi,yi,zi,cmax,cmin,dc,unitLab,titl,pPos,pos)

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
% xh2 = gshhs_c.lon;
% yh2 = gshhs_c.lat;
% lev = gshhs_c.level;
% kl = (lev == 2);
% xh2 = xh2(kl);
% yh2 = yh2(kl);
% xh2 = cell2mat(xh2);
% yh2 = cell2mat(yh2);
% xh2 = xh2(:)';
% yh2 = yh2(:)';

figure(fn)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

% hp2 = patchm(yh2,xh2,-1e+8,[1 1 1]);
% set(hp2,'EdgeColor','none')
%geoshow('landareas.shp','FaceColor',[0.7 0.7 0.7],'EdgeColor','none')


% ------------------------ use instead of scatter -------------------------
%ncolors = 64;
%indCol = dot2color(zi,ncolors,cmax,cmin);
% Col = colormap(linspecer(ncolors));
% % if points do not overlap the following can be used (much faster):
% hold on
% for i=1:ncolors
%     K = indCol == i;
%     plotm(yi(K),xi(K),'o','MarkerSize',1.5,'Color',Col(i,:))
% end
%
% % if points overlap we must use (it will produce same plot as scatterm):
% hold on
%for i=length(zi):-1:1
%     plotm(yi(i),xi(i),'o','MarkerSize',1.5,'Color',Col(indCol(i),:))
% end
% -------------------------------------------------------------------------

% ----------------------------- with scatter ------------------------------
colormap(linspecer);
scatterm(ha,yi,xi,4,zi,'o');
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
cpos(3) = 0.6.*cpos(3);
% -nodisplay mode shifts the colorbar (it is a bug) when printing. Here we
% shift back the colorbar by an amount similar to that done by Matlab  
cpos(2) = cpos(2)+0.765; 
set(hc,'Position',cpos)
set(gca,'position',xpos)
title(titl,'FontSize',10)