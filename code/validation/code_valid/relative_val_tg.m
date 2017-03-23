%relative validation of tide gauges
function relative_val_tg(month)
% input : month: string in format 'mmm-yyyy' ('Mar-2015') representing the
% year up to which the validation is performed

%addpath(genpath('/Volumes/noc/users/cryo/QCV_Cryo2/code/'))
% UPDATED JAN 2017
path1 = '/noc/mpoc/cryo/cryosat/daily_stats/';
path4 = '/noc/mpoc/cryo/cryosat/validation_data_2/'; %gshhs
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';

% path1 = '/scratch/general/cryosat/daily_stats/';
% path4 = '/scratch/general/cryosat/validation_data_2/'; %gshhs
% path_C2Data = '/scratch/general/cryosat/validation_data/';
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/';

% ---------------------- options ------------------------------------------
data_type = 'SIR_GOP_L2';

daten = datenum(month,'mmm-yyyy');
[Y,M] = datevec(daten);
Df = datestr([Y,M,eomday(Y,M),0,0,0],'yyyymmdd');
D0 = '20140415';
detided = 'yes';
smoothed = 'no';
distMax = 60; %maximum distance (km) between altimetry measurement and TG
timeMax = 5; %maximum time difference in days

% fn = dir([dir_out 'tg_data_hourly/OS_UH-*.nc']);
% fn = struct2cell(fn);
% fn = fn(1,:);
% fn = cellfun(@(x) x(10:12),fn,'unif',false);
% nstation = cellfun(@str2num,fn);

fn = dir([path4 'tg_data_hourly/h*.dat']);
fn = struct2cell(fn);
fn = fn(1,:);
fn = cellfun(@(x) x(2:4),fn,'unif',false);
nstation = cellfun(@str2num,fn);

% -------------------------------------------------------------------------

rho = NaN(length(nstation),1);
rho2 = NaN(length(nstation),1); % "2" for series with the tide kept
rms = NaN(length(nstation),1);
rms2 = NaN(length(nstation),1);
nrms = NaN(length(nstation),1);
nrms2 = NaN(length(nstation),1);
lon = NaN(length(nstation),1);
lat = NaN(length(nstation),1);
P = NaN(length(nstation),1);
P2 = NaN(length(nstation),1);
np = NaN(length(nstation),1);
slTg = cell(length(nstation),1);
slTg2 = cell(length(nstation),1);
sla = cell(length(nstation),1);
sla2 = cell(length(nstation),1);
ti = cell(length(nstation),1);
ti2 = cell(length(nstation),1);


for ii = 1:length(nstation)
    disp(ii)
    % ------------------------------- read data -------------------------------
    [xtg,ytg,ztg2,ttg] = read_tgData_ascii(nstation(ii)); % read tide gauge data
    if strcmp(detided,'yes')
        [ztg,~] = detide(ztg2,ttg); % remove tide
    end
    ztg = ztg./1000;
    ztg2 = ztg2./1000;
    
    d0 = datenum(D0,'yyyymmdd');
    df = datenum(Df,'yyyymmdd');
    d = datestr(d0:df,'yyyymmdd');
    ndays = size(d,1);
    files = [repmat(['Stats_' data_type '_'],ndays,1) d repmat('.mat',ndays,1)];
    data = load([path1 files(1,:)]);
    o0 = data.orbits(1);
    data = load([path1 files(end,:)]);
    of = data.orbits(end);
    no = of-o0+1;
    
    v = cell(no,1);
    v2 = cell(no,1);
    vtg = cell(no,1);
    vtg2 = cell(no,1);
    t = cell(no,1);
    x = cell(no,1);
    y = cell(no,1);
    o = cell(no,1);
    kn = 0;
    for i=1:ndays
        dataVal = load([path_C2Data 'dataVal_SIR_GOP_L2_' d(i,:) '.mat']);
        xtemp = dataVal.x_ssha;
        ytemp = dataVal.y_ssha;
        vtemp = dataVal.ssha + dataVal.atm; %tides out, atm in
        ttemp = dataVal.t_ssha;
        otemp = dataVal.orbit_ssha;
        vtemp2 = dataVal.ssh - dataVal.mss; %tides in, atm in
        
        ktemp = Distance(xtg,ytg,xtemp,ytemp)./1000;
        ktemp = find(ktemp < distMax);
        vtgtemp = interp1(ttg,ztg,ttemp);
        vtgtemp2 = interp1(ttg,ztg2,ttemp);
        for j=1:length(ktemp)
            tmin = min(abs(ttemp(j) - ttg));
            if tmin > timeMax
                vtgtemp(j) = NaN;
                vtgtemp2(j) = NaN;
            end
        end
        
        knan = isnan(vtgtemp);
        xtemp(knan) = [];
        ytemp(knan) = [];
        vtemp(knan) = [];
        vtemp2(knan) = [];
        ttemp(knan) = [];
        otemp(knan) = [];
        vtgtemp(knan) = [];
        vtgtemp2(knan) = [];
        
        ou = unique(otemp);
        for j=1:length(ou)
            kn = kn+1;
            ko = find(otemp == ou(j));
            x{kn} = xtemp(ko);
            y{kn} = ytemp(ko);
            t{kn} = ttemp(ko);
            v{kn} = vtemp(ko);
            v2{kn} = vtemp2(ko);
            vtg{kn} = vtgtemp(ko);
            vtg2{kn} = vtgtemp2(ko);
            o{kn} = otemp(ko);
        end
        
    end
    
    
    k = cellfun(@isempty,o);
    x(k) = [];
    y(k) = [];
    v(k) = [];
    v2(k) = [];
    t(k) = [];
    vtg(k) = [];
    vtg2(k) = [];
    
    d = cell(size(x));
    mu = NaN(size(x));
    mu2 = NaN(size(x));
    muTg = NaN(size(x));
    muTg2 = NaN(size(x));
    t2 = NaN(size(x));
    tt2 = NaN(size(x));
    dx = distMax;
    for i=1:length(x)
        d{i} = Distance(x{i},y{i},xtg,ytg)./1000;
        STD = nanstd(v{i}(d{i} < dx));
        STD2 = nanstd(v2{i}(d{i} < dx));
        MU = nanmean(v{i}(d{i} < dx));
        MU2 = nanmean(v2{i}(d{i} < dx));
        vtmp = abs(v{i}-MU);
        vtmp2 = abs(v2{i}-MU2);
        if length(vtmp(d{i} < dx)) > 2
            f1 = vtmp < 3*STD;
        else
            f1 = true(size(vtmp));
        end
        if length(vtmp2(d{i} < dx)) > 2
            f2 = vtmp2 < 3*STD2;
        else
            f2 = true(size(vtmp2));
        end
        
        mu(i) = nanmean(v{i}(d{i} < dx & f1));
        mu2(i) = nanmean(v2{i}(d{i} < dx & f2));
        
        muTg(i) = nanmean(vtg{i}(d{i} < dx & f1));
        muTg2(i) = nanmean(vtg2{i}(d{i} < dx & f2));
        t2(i) = nanmean(t{i}(d{i} < dx & f1));
        tt2(i) = nanmean(t{i}(d{i} < dx & f2));
    end
    
    knan= mu==mu & muTg==muTg;
    np(ii) = sum(knan);
    mu = mu(knan);
    muTg = muTg(knan);
    t2 = t2(knan);
    
    knan= mu2==mu2 & muTg2==muTg2;
    mu2 = mu2(knan);
    muTg2 = muTg2(knan);
    tt2 = tt2(knan);
    
    mu = mu - mean(mu);
    mu2 = mu2 - mean(mu2);
    muTg = muTg -mean(muTg);
    muTg2 = muTg2 -mean(muTg2);
    
    if strcmp(smoothed,'yes')
        yy = smooth(t2,mu,5);
        yy2 = smooth(t2,muTg,5);
        mu = yy;
        muTg = yy2;
    end
    
    disp(length(mu))
    
    lon(ii) = xtg;
    lat(ii) = ytg;
    [r,p] = corrcoef(mu,muTg);
    [r2,p2] = corrcoef(mu2,muTg2);
    rho(ii) = r(1,2);
    rho2(ii) = r2(1,2);
    rms(ii) = std(mu-muTg);
    rms2(ii) = std(mu2-muTg2);
    nrms(ii) = rms(ii)/std(muTg)*100;
    nrms2(ii) = rms2(ii)/std(muTg2)*100;
    P(ii) = p(1,2);
    P2(ii) = p2(1,2);
    slTg{ii} = muTg;
    slTg2{ii} = muTg2;
    sla{ii} = mu;
    sla2{ii} = mu2;
    ti{ii} = t2;
    ti2{ii} = tt2;
    
end

ksignif = P < 0.05;
rms = rms.*100;
ksignif2 = P2 < 0.05;
rms2 = rms2.*100;

mrho = mean(rho);
mrho2 = mean(rho2);
mrms = mean(rms);
mrms2 = mean(rms2);
mnrms = mean(nrms);
mnrms2 = mean(nrms2);



% --------------------- plot correlation map (detided) --------------------
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
scatterm(ha,lat(ksignif),lon(ksignif),50,rho(ksignif),'o','filled');
scatterm(ha,lat(~ksignif),lon(~ksignif),50,'k','o');
textm(-84,60,['mean = ',num2str(mrho,'%.2f')],'FontSize',12)
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
print('-depsc',[dir_out 'Fig_rho_map_tg_detided'],'-r300') % save figure
close all

% -------------------------------------------------------------------------

% ----------------------- plot RMS map (detided) --------------------------
figure(1)
cmin = 0;
cmax = 10;
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
scatterm(ha,lat,lon,50,rms,'o','filled');
textm(-84,60,['mean = ',num2str(mrms,'%.1f'),' cm'],'FontSize',12)
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
set(hc,'YTick',cmin:2:cmax)
set(hc,'FontSize',12)
ylabel(hc,'RMS difference (cm)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_rms_map_tg_detided'],'-r300') % save figure
close all

% -------------------------------------------------------------------------



% ------------------ plot normalized RMS map (detided) --------------------
figure(1)
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
scatterm(ha,lat,lon,50,nrms,'o','filled');
textm(-84,60,['mean = ',num2str(mnrms,'%.1f'),' %'],'FontSize',12)
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
caxis([0 140])
set(hc,'YTick',0:20:140)
set(hc,'FontSize',12)
ylabel(hc,'Normalized RMS difference (%)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_nrms_map_tg_detided'],'-r300') % save figure
close all

% -------------------------------------------------------------------------



% ----------------- plot time series example  1 (detided) -----------------
y1 = slTg{7}.*100;
y2 = sla{7}.*100; % Chatham
t2 = ti{7};

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

hl = legend([hp1 hp2],{'Tide gauge','SSHA CryoSat-2'}, ...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('Detided sea level (cm)','FontSize',10)
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
ylim([ymin-2 ymax+2])
box 
rotateXLabels(gca,45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_ssha_timeSeries_detided_1'],'-r300')
close all

% -------------------------------------------------------------------------



% ----------------- plot time series example 2 (detided) ------------------
y1 = slTg{25}.*100;
y2 = sla{25}.*100; % Corunya
t2 = ti{25};

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

hl = legend([hp1 hp2],{'Tide gauge','SSHA CryoSat-2'},...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('Detided sea level (cm)','FontSize',10)
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
ylim([ymin-2 ymax+2])
box 
rotateXLabels(gca,45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_ssha_timeSeries_detided_2'],'-r300')
close all

% -------------------------------------------------------------------------



% --------------------- plot correlation map (tide in) --------------------
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
scatterm(ha,lat(ksignif2),lon(ksignif2),50,rho2(ksignif2),'o','filled');
scatterm(ha,lat(~ksignif2),lon(~ksignif2),50,'k','o');
textm(-84,60,['mean = ',num2str(mrho2,'%.2f')],'FontSize',12)
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
print('-depsc',[dir_out 'Fig_rho_map_tg'],'-r300') % save figure
close all

% -------------------------------------------------------------------------

% ----------------------- plot RMS map (tide in) --------------------------
figure(1)
cmin = 0;
cmax = 10;
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
scatterm(ha,lat,lon,50,rms2,'o','filled');
textm(-84,60,['mean = ',num2str(mrms2,'%.1f'),' cm'],'FontSize',12)
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
set(hc,'YTick',cmin:2:cmax)
set(hc,'FontSize',12)
ylabel(hc,'RMS difference (cm)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_rms_map_tg'],'-r300') % save figure
close all

% -------------------------------------------------------------------------


% ------------------ plot normalized RMS map (tide in) --------------------
figure(1)
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
scatterm(ha,lat,lon,50,nrms2,'o','filled');
textm(-84,60,['mean = ',num2str(mnrms2,'%.1f'),' %'],'FontSize',12)
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
caxis([0 140])
set(hc,'YTick',0:20:140)
set(hc,'FontSize',12)
ylabel(hc,'Normalized RMS difference (%)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_nrms_map_tg'],'-r300') % save figure
close all

% -------------------------------------------------------------------------


% ----------------- plot time series example 1 (tide in) ------------------
y1 = slTg2{7}.*100;
y2 = sla2{7}.*100; % Chatham
t2 = ti2{7};

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

hl = legend([hp1 hp2],{'Tide gauge','SSHA CryoSat-2'},...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('Total sea level (cm)','FontSize',10)

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
print('-depsc',[dir_out 'Fig_ssha_timeSeries_1'],'-r300')
close all

% -------------------------------------------------------------------------


% ----------------- plot time series example 2 (tide in) ------------------
y1 = slTg2{25}.*100;
y2 = sla2{25}.*100; % Corunya
t2 = ti2{25};

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

hl = legend([hp1 hp2],{'Tide gauge','SSHA CryoSat-2'},...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('Total sea level (cm)','FontSize',10)

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
    'XTickLabel',timeLabel,'YTick',-240:dl:240)
xlim([min(t2)-5 max(t2)+5])
ylim([ymin-5 ymax+5])
box 
rotateXLabels(gca,45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_ssha_timeSeries_2'],'-r300')
close all

% -------------------------------------------------------------------------






