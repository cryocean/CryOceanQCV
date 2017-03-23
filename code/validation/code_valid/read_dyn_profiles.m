
function read_dyn_profiles(month)
% ---------------------------- set paths ----------------------------------
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures

% UPDATED JAN 2017
path1 = '/noc/mpoc/cryo/cryosat/validation_data_2/EN4_TS_profiles/';
path2 = '/noc/mpoc/cryo/cryosat/validation_data_2/';

% path1 = '/scratch/general/cryosat/validation_data_2/EN4_TS_profiles/';
% path2 = '/scratch/general/cryosat/validation_data_2/';
%addpath(genpath('/Users/cryo/MATLAB/'))
zref = 1000;
% -------------------------------------------------------------------------

% load altimetry data
tf = datenum(month,'mmm-yyyy');
[Y,M,~] = datevec(tf);
tf = tf + eomday(Y,M) - 1;
load([path2 'ssha_c2_20140413_5days.mat']);
x = load([path2 'lon_alt.txt']);
y = load([path2 'lat_alt.txt']);
t0 = datenum(2014,04,13);
ta = t0:5:(t0+size(z,1)*5-1); %#ok
ta = ta';
kk = find(ta <= tf);
ta = ta(kk);
z = z(kk,:,:);

load([path2,'gshhs_i.mat']);
x = x(:);
y = y(:);
z = z(:,:);
kland = island(x,y,gshhs_i);
x = x(~kland);
y = y(~kland);
z = z(:,~kland);
k60 = find(abs(y) <= 65);
x = x(k60);
y = y(k60);
z = z(:,k60);

% load gridded steric for climatology
xc = load([path2 'lon_EN4.txt']);
yc = load([path2 'lat_EN4.txt']);
hs = load([path2 'steric_height_1993_2009_EN4_1000m.mat']);
hs = hs.hs;
hs = squeeze(mean(hs(1:204,:,:),1)); % mean steric 1993-2009 (referece period CryoSat)
xc(xc>180) = xc(xc>180)-360;
[xc,ix] = sort(xc);
hs = hs(ix,:);
[yc,xc] = meshgrid(yc,xc);
knanc = find(~isnan(hs));
xc = xc(knanc);
yc = yc(knanc);
hs = hs(knanc);
Fc = TriScatteredInterp(xc,yc,hs);

% save profile data in a map container (each key is a profile id)
fn = dir([path1 'dataArgo*' num2str(zref) '.mat']);
fn = fn(4:end);
map = containers.Map('KeyType','char', 'ValueType','any');
for i=1:length(fn)
    disp(i)
    datai = load([path1 fn(i).name]);
    datai = datai.dataArgo;
    idi = {datai(:).id};
    for j=1:length(datai);
        xp = datai(j).lon;
        yp = datai(j).lat;
        tp = datai(j).time;
        hp = datai(j).steric;
        for k=1:length(xp)
            if isKey(map,idi(j))
                tmp = [xp(k); yp(k); tp(k); hp(k)];
                map(idi{j}) = [map(idi{j}) tmp];
            else
                map(idi{j}) = [xp(k); yp(k); tp(k); hp(k)];
            end
        end
        
    end
    
end

key = keys(map);
Hp = cell(length(key),1);
Ha = cell(length(key),1);
Hclim = cell(length(key),1);
Xp = cell(length(key),1);
Yp = cell(length(key),1);
Tp = cell(length(key),1);

Fs = cell(length(ta),1);
for i=1:length(ta)
    Fs{i} = TriScatteredInterp(x,y,z(i,:)','linear'); %spatial interpolation
end
for i=1:length(key)
    disp(i)
    prof = map(key{i});
    xp = prof(1,:);
    yp = prof(2,:);
    tp = prof(3,:);
    hp = prof(4,:);
    ha = zeros(1,length(xp));
    hmean = zeros(1,length(xp));
    if abs(yp(1)) < 60
        for j=1:length(xp)
            hmean(j) = Fc(xp(j),yp(j)); % climatology interpolated to profile 
            htmp = zeros(size(z,1),1);
            for k=1:size(z,1)
                htmp(k) = Fs{k}(xp(j),yp(j)); %spatial interpolation
            end
            ha(j) = interp1(ta,htmp,tp(j)); % temporal interpolation
        end
        try
            k0 = find(~isnan(ha),1,'first');
            kf = find(~isnan(ha),1,'last');
            Hp{i} = hp(k0:kf);
            Ha{i} = ha(k0:kf);
            Hclim{i} = hmean(k0:kf);
            Xp{i} = xp(k0:kf);
            Yp{i} = yp(k0:kf);
            Tp{i} = tp(k0:kf);
            
        catch %#ok
            continue
        end
    end
end

K = cellfun(@(x,y) sum(~isnan(x) & ~isnan(y)) > 2,Hp,Ha); %profiles with at least 2 obs
zp = Hp(K);
za = Ha(K);
zc = Hclim(K);
tp = Tp(K);
lon = cellfun(@mean,Xp(K));
lat = cellfun(@mean,Yp(K));
xp = Xp(K);
yp = Yp(K);

[rho,p] = cellfun(@(x,y) nancorrcoef(x,y),zp,za);
nrms = cellfun(@(x,y,z) sqrt(nanmean((x-y-z).^2))/(max(y)-min(y))*100,zp,za,zc);
bias = cellfun(@(x,y,z) nanmean(y-(x-z)),zp,za,zc);
load([path2, 'gshhs_c.mat'])
kLand = island(lon(:),lat(:),gshhs_c); %profiles over land
Ks = find(rho > 0 & p<0.05 & ~kLand');


% ---------------------- plot correlation map -----------------------------
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

figure(1)
cmin = 0;
cmax = 1;
pPos = [13.9 9.4];
pos = [10.5 9.2];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

colormap(linspecer);
scatterm(ha,lat(Ks),lon(Ks),30,rho(Ks),'o','filled');
framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
tightmap
%title('Location of tide gauge stations used in the validation','FontSize',12)

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([cmin cmax])
set(hc,'YTick',cmin:0.2:cmax)
set(hc,'FontSize',12)
ylabel(hc,'Correlation') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_rho_map_dyn'],'-r300') % save figure
close all
% -------------------------------------------------------------------------


% ----------------------- plot normalized RMSD ----------------------------
figure(1)
cmin = 0;
cmax = 100;
pPos = [13.9 9.4];
pos = [10.5 9.2];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

colormap(linspecer);
scatterm(ha,lat(Ks),lon(Ks),30,nrms(Ks),'o','filled');
framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
tightmap
%title('Location of tide gauge stations used in the validation','FontSize',12)

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([cmin cmax])
set(hc,'YTick',cmin:20:cmax)
set(hc,'FontSize',12)
ylabel(hc,'Normalized RMSD (%)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)
mu = nanmean(nrms(Ks));
textm(-85,60,['Mean = ',sprintf('%3.1f',mu),'%'],'FontSize',12)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_NRMSD_map_dyn'],'-r300') % save figure
close all
% -------------------------------------------------------------------------



% ----------------------- plot bias (SLA - DYN) ---------------------------
figure(1)
cmin = -15;
cmax = 5;
pPos = [13.9 9.4];
pos = [10.5 9.2];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

colormap(linspecer);
scatterm(ha,lat(Ks),lon(Ks),30,bias(Ks).*100,'o','filled');
framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
tightmap
%title('Location of tide gauge stations used in the validation','FontSize',12)

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([cmin cmax])
set(hc,'YTick',cmin:5:cmax)
set(hc,'FontSize',12)
ylabel(hc,'Bias: SSHA - Steric height (cm)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)
mu = nanmean(bias(Ks)).*100;
textm(-85,60,['Mean = ',sprintf('%3.1f',mu),' cm'],'FontSize',12)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_bias_map_dyn'],'-r300') % save figure
close all
% -------------------------------------------------------------------------


% --------------------- plot time series examples -------------------------
% pick 2 time series randomly among those with rho>0.65 and nrms<30
K2 = find(rho(Ks)>0.65 & nrms(Ks)<30);
ik = datasample(1:length(K2),2,'Replace',false);

% first profile
y1 = (zp{Ks(K2(ik(1)))}-nanmean(zp{Ks(K2(ik(1)))})).*100;
y2 = (za{Ks(K2(ik(1)))}-nanmean(za{Ks(K2(ik(1)))})).*100;
t2 = tp{Ks(K2(ik(1)))};
lon1 = xp{Ks(K2(ik(1)))};
lat1 = yp{Ks(K2(ik(1)))};
idProfile = key(Ks(K2(ik(1))));
% knan = isnan(y1) | isnan(y2);
% y1(knan) = NaN;
% y2(knan) = NaN;

pPos = [11.5 6.7];
pos = [8.7 5];
figure
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.5 1.5 pos])
rgb1 = [215/255 75/255 75/255];
rgb2 = [120/255 121/255 116/255];
hold on
hp1 = plot(t2,y1,'-o','color',rgb1,'LineWidth',0.7,'MarkerSize',2, ...
    'MarkerFaceColor',rgb1);
hp2 = plot(t2,y2,'-o','color',rgb2,'LineWidth',0.7,'MarkerSize',2, ...
    'MarkerFaceColor',rgb2);

hl = legend([hp1 hp2],{'Argo steric','SSHA CryoSat-2'},...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('SSHA/Steric (cm)','FontSize',10)
%dx = ceil(length(t2)/5);
dx = round((t2(end)-t2(1))/8);
dtt = [t2(1):dx:t2(end),t2(end)];
if dtt(end)-dtt(end-1) < dx/2.5 || dtt(end) == dtt(end-1)
    dtt(end-1) = [];
end

ymin = min([y1(:);y2(:)]);
ymax = max([y1(:);y2(:)]);
dtmp = ymax - ymin + 4;
dl = 0;
c = 1;
while dl == 0
    dl = round(dtmp*c/5)/c;
    c = c*10;
end

timeLabel = cellstr(datestr(dtt,'dd/mm/yy'));
set(gca,'LineWidth',1,'FontSize',10,'XTick',dtt, ...
    'XTickLabel',timeLabel,'YTick',-200:dl:200)
xlim([min(t2)-5 max(t2)+5])
ylim([ymin-5 ymax+5])
box 
rotateXLabels(gca,45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_ssha_steric_timeSeries_',num2str(1)],'-r300')
close all


% plot small map showing location of profiles
figure(1)
pPos = [13.5 9.5];
pos = [10 8];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 0.8 pos])
ha = axesm('MapProjection','gstereo','MapLatLimit',[median(lat1)-7 median(lat1)+7],...
    'MapLonLimit',[median(lon1)-10 median(lon1)+10]);
hp = patchm(yh,xh,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

colormap(linspecer);
t0 = t2(1);
dt = t2-t0;
scatterm(lat1,lon1,30,dt,'d','filled')

textm(median(lat1)+6,median(lon1)-9,['Argo: ' idProfile],'FontSize',12)

framem
mlabel('MLabelParallel','south','MLabelLocation',-180:5:180)
plabel('PLabelLocation',-90:4:90)
setm(ha,'FontSize',12)
tightmap
title('Location of profile over time','FontSize',12)

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([0 max(dt)])
dtmp = max(dt);
dy = 0;
ctmp = 1;
while dy == 0
    dy = round(dtmp*ctmp/5)/ctmp;
    ctmp = ctmp*10;
end
set(hc,'YTick',0:dy:max(dt))
set(hc,'FontSize',12)
ylabel(hc,'Days since first sample') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.711;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_profile_location_',num2str(1)],'-r300') % save figure
close all




% second profile
y1 = (zp{Ks(K2(ik(2)))}-nanmean(zp{Ks(K2(ik(2)))})).*100;
y2 = (za{Ks(K2(ik(2)))}-nanmean(za{Ks(K2(ik(2)))})).*100;
t2 = tp{Ks(K2(ik(2)))};
lon1 = xp{Ks(K2(ik(2)))};
lat1 = yp{Ks(K2(ik(2)))};
idProfile = key(Ks(K2(ik(2))));
knan = isnan(y1) | isnan(y2);
y1(knan) = NaN;
y2(knan) = NaN;

pPos = [11.5 6.7];
pos = [8.7 5];
figure
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.5 1.5 pos])
rgb1 = [215/255 75/255 75/255];
rgb2 = [120/255 121/255 116/255];
hold on
hp1 = plot(t2,y1,'-o','color',rgb1,'LineWidth',0.7,'MarkerSize',2, ...
    'MarkerFaceColor',rgb1);
hp2 = plot(t2,y2,'-o','color',rgb2,'LineWidth',0.7,'MarkerSize',2, ...
    'MarkerFaceColor',rgb2);

hl = legend([hp1 hp2],{'Argo steric','SSHA CryoSat-2'},...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('SSHA/Steric (cm)','FontSize',10)
%dx = ceil(length(t2)/5);
dx = round((t2(end)-t2(1))/8);
dtt = [t2(1):dx:t2(end),t2(end)];
if dtt(end)-dtt(end-1) < dx/2.5 || dtt(end) == dtt(end-1)
    dtt(end-1) = [];
end

ymin = min([y1(:);y2(:)]);
ymax = max([y1(:);y2(:)]);
dtmp = ymax - ymin + 4;
dl = 0;
c = 1;
while dl == 0
    dl = round(dtmp*c/5)/c;
    c = c*10;
end

timeLabel = cellstr(datestr(dtt,'dd/mm/yy'));
set(gca,'LineWidth',1,'FontSize',10,'XTick',dtt, ...
    'XTickLabel',timeLabel,'YTick',-200:dl:200)
xlim([min(t2)-5 max(t2)+5])
ylim([ymin-5 ymax+5])
box 
rotateXLabels(gca,45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_ssha_steric_timeSeries_',num2str(2)],'-r300')
close all


% plot small map showing location of profiles
figure(1)
pPos = [13.5 9.5];
pos = [10 8];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 0.8 pos])
ha = axesm('MapProjection','gstereo','MapLatLimit',[median(lat1)-7 median(lat1)+7],...
    'MapLonLimit',[median(lon1)-10 median(lon1)+10]);
hp = patchm(yh,xh,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

colormap(linspecer);
t0 = t2(1);
dt = t2-t0;
scatterm(lat1,lon1,30,dt,'d','filled')

textm(median(lat1)+6,median(lon1)-9,['Argo: ' idProfile],'FontSize',12)

framem
mlabel('MLabelParallel','south','MLabelLocation',-180:5:180)
plabel('PLabelLocation',-90:4:90)
setm(ha,'FontSize',12)
tightmap
title('Location of profile over time','FontSize',12)

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([0 max(dt)])
dtmp = max(dt);
dy = 0;
ctmp = 1;
while dy == 0
    dy = round(dtmp*ctmp/5)/ctmp;
    ctmp = ctmp*10;
end
set(hc,'YTick',0:dy:max(dt))
set(hc,'FontSize',12)
ylabel(hc,'Days since first sample') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.711;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_profile_location_',num2str(2)],'-r300') % save figure
close all
% -------------------------------------------------------------------------

