function plot_mode_mask(seas,day_date)
dir1 = '/noc/users/cryo/QCV_Cryo2/code/mode_mask/';
path4 = '/noc/users/cryo/QCV_Cryo2/code/';
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/';

if datenum(day_date,'yyyymmdd') >= datenum('20160307','yyyymmdd')
    mask_v = '8';
elseif datenum(day_date,'yyyymmdd') >= datenum('20151214','yyyymmdd')
    mask_v = '7';
elseif datenum(day_date,'yyyymmdd') >= datenum('20141006','yyyymmdd')
    mask_v = '6';
elseif datenum(day_date,'yyyymmdd') >= datenum('20111231','yyyymmdd')
    mask_v = '5';
else
    mask_v = '2';
end

load([dir1 'seasMaskAllModes_3_' mask_v '.mat'])
mask = seasMaskModes.(['seas_' sprintf('%1.2d',seas)]);
Mode = mask.mode;
X = mask.X;
Y = mask.Y;
p = mask.priority;
ksar = cellfun(@strcmp,Mode,repmat({'SAR'},1,length(X)));
ksin = cellfun(@strcmp,Mode,repmat({'SIN'},1,length(X)));
klrm = cellfun(@strcmp,Mode,repmat({'LRM'},1,length(X)));
ksar1 = p <= 6 & ksar;
ksar2 = p > 6 & ksar;
ksin = find(ksin);
ksar1 = find(ksar1);
ksar2 = find(ksar2);
klrm = find(klrm);

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

pPos = [16.8 11.04];
pos = [12.6 10.8];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

hold on
for i=1:length(ksar1)
    h1 = patchm(Y{ksar1(i)},X{ksar1(i)},'Facecolor',[1 197/255 108/255], ...
        'EdgeColor','none','LineStyle','none'); %#ok
end

hold on
for i=1:length(ksin)
    h2 = patchm(Y{ksin(i)},X{ksin(i)},'Facecolor',[110/255 197/255 233/255], ...
        'EdgeColor','none','LineStyle','none'); %#ok
end

hold on
for i=1:length(ksar2)
    h1 = patchm(Y{ksar2(i)},X{ksar2(i)},'Facecolor',[1 197/255 108/255], ...
        'EdgeColor','none','LineStyle','none'); %#ok
end

hold on
for i=1:length(klrm)
    h3 = patchm(Y{klrm(i)},X{klrm(i)},'Facecolor',[1 89/255 89/255], ...
        'EdgeColor','none','LineStyle','none'); %#ok
end
framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
geoshow('landareas.shp','FaceColor','none','EdgeColor','black')
tightmap
%legend([h1 h2 h3],'SAR','SARIN','LRM')
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_modeMask'],'-r300') % save figure
close all