function cwsp_val_v2(month)
% validation of CryoSat-2 wind speed against NDBC buoys (continuous wsp)
% NDBC data: http://www.ndbc.noaa.gov/measdes.shtml
%            http://www.ndbc.noaa.gov/data/
% time difference is 5 minutes
% NDBC provides data every 10 minutes
%
% input : month: string in format 'mmm-yyyy' ('Mar-2015') representing the
% year up to which the validation is performed

% UPDATED HERE JAN 2017
path2Code = '/noc/mpoc/cryo/cryosat/validation_data_2/buoys/';
path2 = '/noc/mpoc/cryo/cryosat/validation_data_2/'; %gshhs
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path3 = '/noc/mpoc/cryo/cryosat/daily_stats/';

%         path_C2Data = '/scratch/general/cryosat/validation_data/';
%         path3 = '/scratch/general/cryosat/daily_stats/';


% path2Code = '/scratch/general/cryosat/validation_data_2/buoys/';
% path2 = '/scratch/general/cryosat/validation_data_2/'; %gshhs

dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures
pathOut = [path2Code,'c2BuoyWsp/'];
minDist = 20; % minimum distance from coast
maxDist = 10; % maximum distance between buoy and satellite measurement
data_type = 'SIR_GOP_L2';
date0 = datenum('20140412','yyyymmdd');
datef = datenum(month,'mmm-yyyy');
days = date0:datef;
days = cellstr(datestr(days,'yyyymm'));
days = unique(days);


Lon = cell(length(days),1);
Lat = cell(length(days),1);
Month = cell(length(days),1);
V_c2 = cell(length(days),1);
V_buoy = cell(length(days),1);
for mm=1:length(days)
    disp(mm)
    if exist([pathOut,'c2BuoyData','_',days{mm},'_',num2str(minDist),'km.mat'], 'file') == 2
        load([pathOut,'c2BuoyData','_',days{mm},'_',num2str(minDist),'km.mat']);
    else
        c2BuoyData = struct('lon',[],'lat',[],'month',[],'v_c2',[],'v_buoy',[]);
        date0 = datenum(days{mm},'yyyymm');
        [Y,M] = datevec(date0);
        datef = date0+eomday(Y,M)-1;
        if Y == 2014 && M == 4;
            date0 = datenum(2014,4,12); % no CryoSat2 data prior to that date
        end
        D0 = datestr(date0,'yyyymmdd');
        Df = datestr(datef,'yyyymmdd');
        
        
        % read station id and coordinates for all stations in the data base
        fidId = fopen([path2Code 'stationId.txt']);
        stationId = textscan(fidId,'%s');
        stationId = stationId{1};
        stationId = lower(stationId); %to lowercase
        lon = load([path2Code 'lonStation.txt']);
        lat = load([path2Code 'latStation.txt']);
        fclose(fidId);
        
        if Y == 2014
            path1 = [path2Code num2str(Y) '_cwsp/'];
        else
            path1 = [path2Code num2str(Y) '_cwsp/',datestr(date0,'mmm'),'/'];
        end
        
        % retrieve id of buoys with swh and wsp data
        fn = dir([path1 '*.txt']);
        fn = struct2cell(fn);
        fn = fn(1,:)';
        idData = regexp(fn,'\.','split');
        idData = cellfun(@(x) x(1),idData); %id's of stations with data
        
        % select stations that have data
        [~,kid,~] = intersect(stationId,idData);
        %stationId = stationId(kid);
        lon = lon(kid);
        lat = lat(kid);
        stationId = stationId(kid);
        [~,kid] = setdiff(idData,stationId);
        fn(kid) = [];
        idData(kid) = [];
        [idData,kid] = sort(idData);
        fn = fn(kid);
        [stationId,kid] = sort(stationId);
        lon = lon(kid);
        lat = lat(kid);
        
        % read data
        % read data
        wsp = cell(length(fn),1);
        t = cell(length(fn),1);
        format1 = repmat('%d',[1 5]);
        format2 = repmat('%f',[1 5]);
        for i=1:length(fn)
            fid = fopen([path1 fn{i}],'r'); %open data file
            header = fgetl(fid);
            header = regexp(header,'\s+','split');
            dum = fgetl(fid); %#ok discard 2nd row of units
            [~,k,~] = intersect(header,{'WSPD'}); %identify column
            data = fscanf(fid,[format1 format2]) ; %read data
            data = reshape(data,10,length(data)/10);
            tmp1 = data(k(1),:);
            tmp1(tmp1 > 95) = NaN;
            wsp{i} = tmp1;
            t{i} = ([data(1:5,:); zeros(1,size(data,2))])'; %add zeros for seconds
            fclose(fid); %close file
        end
        
        kwsp = cellfun(@(x) any(~isnan(x)),wsp);
        
        wsp = wsp(kwsp);
        twsp = t(kwsp);
        lonwsp = lon(kwsp);
        latwsp = lat(kwsp);
        
        twsp = cellfun(@datenum,twsp,'unif',0);
        
        % prepare data to compute distance from coast
        load([path2 'gshhs_c.mat'])
        
        xh = gshhs_c.lon;
        yh = gshhs_c.lat;
        lev = gshhs_c.level;
        
        %boundary between Antarctica grounding-line and ocean (level 6) rather than
        %Antarctica ice and ocean (level 5)
        kl = (lev == 1 | lev == 2 | lev == 6); %include lakes and enclosed seas
        xh = xh(kl);
        yh = yh(kl);
        xh = cell2mat(xh);
        yh = cell2mat(yh);
        xh = xh(:);
        yh = yh(:);
        xh = xh(~isnan(xh));
        yh = yh(~isnan(yh));
        [yh,IX] = sort(yh);
        xh = xh(IX);
        xh(xh > 180) = xh(xh > 180)-360;
        
        % reject stations that are closer than 20km from coast
        dwsp = distance2Coast(lonwsp(:),latwsp(:),xh(:),yh(:),1000)./1000;
        kwsp = find(dwsp > minDist);
        wsp = wsp(kwsp);
        twsp = twsp(kwsp);
        lonwsp = lonwsp(kwsp);
        latwsp = latwsp(kwsp);
        
        
        nstation = length(twsp);
        %         path_C2Data = '/scratch/general/cryosat/validation_data/';
        %         path3 = '/scratch/general/cryosat/daily_stats/';
        
        rms = NaN(nstation,1);
        v_c2 = cell(nstation,1);
        t_c2 = cell(nstation,1);
        v_buoy = cell(nstation,1);
        lon = NaN(nstation,1);
        lat = NaN(nstation,1);
        np = NaN(nstation,1);
        disp(nstation)
        for ii = 1:nstation
            disp(ii)
            % ------------------------------- read data ---------------------------
            ztg = wsp{ii};
            ttg = twsp{ii};
            xtg = lonwsp(ii);
            ytg = latwsp(ii);
            
            if length(ztg) > 2
                
                d0 = datenum(D0,'yyyymmdd');
                df = datenum(Df,'yyyymmdd');
                d = datestr(d0:df,'yyyymmdd');
                %dates = d0:df;
                ndays = size(d,1);
                files = [repmat(['Stats_' data_type '_'],ndays,1) d repmat('.mat',ndays,1)];
                data = load([path3 files(1,:)]);
                o0 = data.orbits(1);
                data = load([path3 files(end,:)]);
                of = data.orbits(end);
                no = of-o0+1;
                
                v = cell(no,1);
                vtg = cell(no,1);
                t = cell(no,1);
                x = cell(no,1);
                y = cell(no,1);
                o = cell(no,1);
                kn = 0;
                for i=1:ndays
                    
                    dataVal = load([path_C2Data 'dataVal_SIR_GOP_L2_' d(i,:) '.mat']);
                    xtemp = dataVal.x_wsp;
                    ytemp = dataVal.y_wsp;
                    vtemp = dataVal.wsp;
                    ttemp = dataVal.t_wsp;
                    otemp = dataVal.orbit_wsp;
                    
                    ktemp = Distance(xtg,ytg,xtemp,ytemp)./1000;
                    ktemp = find(ktemp < maxDist);
                    xtemp = xtemp(ktemp);
                    ytemp = ytemp(ktemp);
                    vtemp = vtemp(ktemp);
                    ttemp = ttemp(ktemp);
                    otemp = otemp(ktemp);
                    [~,ix] = unique(ttg,'last');
                    vtgtemp = interp1(ttg(ix),ztg(ix),ttemp);
                    dtime = zeros(length(ktemp),1);
                    for j=1:length(ktemp)
                        dtime(j) = min(abs(ttemp(j) - ttg));
                    end
                    
                    kt = find(dtime < 1/(12*24)); %5 minutes time difference
                    if isempty(kt)
                        continue
                    end
                    
                    xtemp = xtemp(kt);
                    ytemp = ytemp(kt);
                    ttemp = ttemp(kt);
                    vtemp = vtemp(kt);
                    otemp = otemp(kt);
                    vtgtemp = vtgtemp(kt);
                    
                    ou = unique(otemp); %loop over each pass
                    for j=1:length(ou)
                        kn = kn+1;
                        ko = find(otemp == ou(j));
                        x{kn} = xtemp(ko);
                        y{kn} = ytemp(ko);
                        t{kn} = ttemp(ko);
                        v{kn} = vtemp(ko);
                        vtg{kn} = vtgtemp(ko);
                        o{kn} = otemp(ko);
                    end
                    
                end
                
                k = cellfun(@isempty,o);
                x(k) = [];%#ok
                y(k) = [];%#ok
                v(k) = [];
                t(k) = [];
                o(k) = [];%#ok
                vtg(k) = [];
                
                mu = cellfun(@nanmean,v);
                mu2 = cellfun(@nanmean,vtg);
                t2 = cellfun(@nanmean,t);
                
                knan= mu==mu & mu2==mu2;
                np(ii) = sum(knan);
                mu = mu(knan);
                mu2 = mu2(knan);
                t2 = t2(knan);
                
                rms(ii) = std(mu-mu2);
                lon(ii) = xtg;
                lat(ii) = ytg;
                v_c2{ii} = mu;
                t_c2{ii} = t2;
                v_buoy{ii} = mu2;
                
            end
            
        end
        
        kempty = cellfun(@isempty,v_c2);
        
        v_c2(kempty) = [];
        v_buoy(kempty) = [];
        lon(kempty) = [];
        lat(kempty) = [];
        t_c2(kempty) = [];%#ok
        
        c2BuoyData.lon = lon;
        c2BuoyData.lat = lat;
        c2BuoyData.month = days{mm};
        c2BuoyData.v_c2 = v_c2;
        c2BuoyData.v_buoy = v_buoy;%#ok
        
        save([pathOut,'c2BuoyData','_',days{mm},'_',num2str(minDist),'km.mat'],...
            '-struct','c2BuoyData');
    end
    
    Lon{mm} = lon;
    Lat{mm} = lat;
    Month{mm} = days{mm};
    V_c2{mm} = v_c2;
    V_buoy{mm} = v_buoy;
    
end


% plot results
v1 = V_c2{end};
v2 = V_buoy{end};
v1 = cell2mat(v1);
v2 = cell2mat(v2);

kzero = v2 == 0;
v1(kzero) = [];
v2(kzero) = [];

lon = [];
lat = [];
month = [];
tmp = V_c2{end};
n = cellfun(@length,tmp);
for j=1:length(n)
    lon = [lon;repmat(Lon{end}(j),n(j),1)];%#ok
    lat = [lat;repmat(Lat{end}(j),n(j),1)];%#ok
    month = [month;repmat(Month{end},n(j),1)]; %#ok
end
lon = lon(~kzero);
lat = lat(~kzero);

% reject buoys located in the Great Lakes
K = find(~(lon>-93 & lon<-75 & lat>40 & lat<50));
lon = lon(K);
lat = lat(K);
v1 = v1(K);
v2 = v2(K);

% v1 = cellfun(@nanmean,v_c2);
% v2 = cellfun(@nanmean,v_buoy);

% calculate OLS
XYLIM = [0 25];
tt = (XYLIM(1):0.1:XYLIM(2))';
X = [ones(size(tt)) tt];
stats = regstats(v1,v2,'linear');
rols = X*stats.beta;
bols = stats.beta;
eols = stats.tstat.se;
R2 = stats.rsquare;
rmsols = sqrt(stats.mse);

% calculate Robustfit
[brob,stats] = robustfit(v2,v1);
rrob = X*brob;
erob = stats.se;
rmsrob = stats.s;

RMS = zeros(length(Lon),1);
Dif = zeros(length(Lon),1);
TAU = zeros(length(Lon),1);
for i=1:length(Lon)
    vtmp1 = V_c2{i};
    n = cellfun(@length,vtmp1);
    vtmp2 = V_buoy{i};
    vtmp1 = cell2mat(vtmp1);
    vtmp2 = cell2mat(vtmp2);
    
    lon = [];
    lat = [];
    for j=1:length(n)
        lon = [lon;repmat(Lon{i}(j),n(j),1)];%#ok
        lat = [lat;repmat(Lat{i}(j),n(j),1)];%#ok
    end
    
    k = vtmp1~=0 & vtmp2 ~=0 & ~(lon>-93 & lon<-75 & lat>40 & lat<50);
    [brob2,stats2] = robustfit(vtmp2(k),vtmp1(k));
    RMS(i) = stats2.s;
    Dif(i) = mean(vtmp1 - vtmp2);
    TAU(i) = brob2(2);
end

K = find(~(lon>-93 & lon<-75 & lat>40 & lat<50));
lon = lon(K);
lat = lat(K);


% ------------------------plot RMS over time ------------------------------
Month2 = Month;
Month2 = datenum(Month2,'yyyymm');
Month2 = datestr(Month2,'mmm-yyyy');
Month2 = cellstr(Month2);
if(length(Month2)>11)
    dx = round(length(Month2)/6);
else
    dx = 2;
end
if length(Month2) > 7
    xti = [1:dx:(length(Month2)-2) length(Month2)];
    Month2 = {Month2{1:dx:end-2},Month2{end}};
else
    xti = 1:length(Month2);
end

ymax = max(RMS) + (max(RMS) - min(RMS))/20;
ymin = min(RMS) - (max(RMS) - min(RMS))/20;
dtmp = ymax - ymin;
dy = 0;
c = 1;
while dy == 0
    dy = round(dtmp*c/5)/c;
    c = c*10;
end

% plot bias
figure(1)
pPos = [15 15];
pos = [12 6];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.7 2.6 pos])

hold on
plot(RMS,'-o','LineWidth',1,'color','k','MarkerFaceColor','k','MarkerSize',4)
plot(length(RMS),RMS(end),'o','LineWidth',1,'color','k',...
    'color','r','MarkerFaceColor','r','MarkerSize',4)
box on
xlim([0 length(RMS)+1])
ylabel('RMS wind speed (m/s)')
set(gca,'FontSize',12,'XTick',xti,'XTickLabel',Month2,'YTick',0:dy:5)
ylim([ymin ymax])
rotateXLabels(gca(),45)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_cwsp_RMS'],'-r300') % save figure
close all
% -------------------------------------------------------------------------

% ------------------------plot slope over time ------------------------------
ymax = max(TAU) + (max(TAU) - min(TAU))/20;
ymin = min(TAU) - (max(TAU) - min(TAU))/20;
dtmp = ymax - ymin;
dy = 0;
c = 1;
while dy == 0
    dy = round(dtmp*c/5)/c;
    c = c*10;
end

% plot bias
figure(1)
pPos = [15 15];
pos = [12 6];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.7 2.6 pos])

hold on
plot(TAU,'-o','LineWidth',1,'color','k','MarkerFaceColor','k','MarkerSize',4)
plot(length(TAU),TAU(end),'o','LineWidth',1,'color','k',...
    'color','r','MarkerFaceColor','r','MarkerSize',4)
box on
xlim([0 length(TAU)+1])
ylabel('Slope wind Speed')
set(gca,'FontSize',12,'XTick',xti,'XTickLabel',Month2,'YTick',0:dy:2)
ylim([ymin ymax])
rotateXLabels(gca(),45)
%ylim([0.9 1.5])
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_cwsp_TAU'],'-r300') % save figure
close all
% -------------------------------------------------------------------------


%plot for OLS approach
R2 = sprintf('%4.2f',R2);
rmsols = sprintf('%4.2f',rmsols);
sols = sprintf('%4.2f',bols(2));
iols = sprintf('%4.2f',bols(1));
esols = sprintf('%4.2f',eols(2));
eiols = sprintf('%4.2f',eols(1));

rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
pPos = [10 10];
pos = [8 8];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 1 pos])
scatter(v2,v1,15,rgb1,'filled')
ylabel(gca,'CryoSat-2 wind speed (m/s)','FontSize',10)
xlabel(gca,'Buoy wind speed (m/s)','FontSize',10)
set(gca,'LineWidth',1,'FontSize',10,'XTick',0:5:25,'YTick',0:5:25)
xlim([0 25])
ylim([0 25])
hold on
plot(tt,rols,'-','color',rgb2,'LineWidth',1.5)
plot(-1e2:1e2,-1e2:1e2,'--','color','black','LineWidth',1)
text(12,6,['R{^2} = ' R2])
text(12,4.5,['rms = ' rmsols ' m/s'])
text(12,3.0,['slope = ' sols ' \pm ' esols])
text(12,1.5,['intcp = ' iols ' \pm ' eiols ' m/s'])
box on
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_cwsp_ols'],'-r300') % save figure
close all


%plot for robustfit approach
rmsrob = sprintf('%4.2f',rmsrob);
srob = sprintf('%4.2f',brob(2));
irob = sprintf('%4.2f',brob(1));
esrob = sprintf('%4.2f',erob(2));
eirob = sprintf('%4.2f',erob(1));

rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
pPos = [10 10];
pos = [8 8];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.2 1 pos])
scatter(v2,v1,15,rgb1,'filled')
ylabel(gca,'CryoSat-2 wind speed (m/s)','FontSize',10)
xlabel(gca,'Buoy wind speed (m/s)','FontSize',10)
set(gca,'LineWidth',1,'FontSize',10,'XTick',0:5:25,'YTick',0:5:25)
xlim([0 25])
ylim([0 25])
hold on
plot(tt,rrob,'-','color',rgb2,'LineWidth',1.5)
plot(-1e2:1e2,-1e2:1e2,'--','color','black','LineWidth',1)
text(12,6,['rms = ' rmsrob ' m/s'])
text(12,4.5,['slope = ' srob ' \pm ' esrob])
text(12,3.0,['intcp = ' irob ' \pm ' eirob ' m/s'])
box on
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_cwsp_robust'],'-r300') % save figure
close all

% plot differences
units = ' m/s';
pPos = [12 11];
pos = [8 8];
d = v1-v2;
dm = mean(d);
DM = sprintf('%4.2f',dm);
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.4 1 pos])
plot(d,'x','color',rgb1)
ylabel(gca,['CryoSat-2 - buoys (' units(2:end) ')'],'FontSize',10)
xlabel(gca,'Measurement number','FontSize',10)
xx = 0:(length(v1)+1);
xlim([0 length(v1)+1])
hold on
plot(xx,dm*ones(1,length(xx)),'-','color',rgb2,'LineWidth',1.5)
plot(xx,zeros(1,length(xx)),'--','color','black','LineWidth',1.5)

YLIM = get(gca,'Ylim');
XLIM = get(gca,'XLIM');
ypos = (YLIM(2)-YLIM(1))*0.95+YLIM(1);
xpos = (XLIM(2)-XLIM(1))*0.05+XLIM(1);
text(xpos,ypos,['mean = ' DM units])
set(gca,'LineWidth',1,'FontSize',10)
box on
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_cwsp_diff'],'-r300') % save figure
close all


% ---------------------------- plot map location --------------------------
load([path2 'gshhs_c.mat'])
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
pPos = [13.9 9.4];
pos = [10.5 9.2];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])


scatterm(ha,lat,lon,15,'k','o','filled');
framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
tightmap
title('Location of the buoys used in the validation','FontSize',12)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_buoy_cwsp_location'],'-r300') % save figure
close all
% -------------------------------------------------------------------------

