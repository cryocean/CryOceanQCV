function swh_wsp_val_v2(month)
% validation of CryoSat-2 SWH and wind speed against NDBC buoys
% NDBC data: http://www.ndbc.noaa.gov/measdes.shtml
%            http://www.ndbc.noaa.gov/data/
% time difference is 30 minutes
% NDBC provides hourly data
%
% input : month: string in format 'mmm-yyyy' ('Mar-2015') representing the
% year up to which the validation is performed
% UPDATED JAN 2017
path2Code = '/noc/mpoc/cryo/cryosat/validation_data_2/buoys/';
path2 = '/noc/mpoc/cryo/cryosat/validation_data_2/'; %gshhs
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path3 = '/noc/mpoc/cryo/cryosat/daily_stats/'; %CryoSat-2 data
% path2Code = '/scratch/general/cryosat/validation_data_2/buoys/';
% path2 = '/scratch/general/cryosat/validation_data_2/'; %gshhs
% path_C2Data = '/scratch/general/cryosat/validation_data/';
% path3 = '/scratch/general/cryosat/daily_stats/'; %CryoSat-2 data



dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures
pathOut = [path2Code,'c2BuoySwh/'];
minDist = 20; % minimum distance from coast
maxDist = 20; % maximum distance between buoy and satellite measurement
vari = 'swh'; % "swh" for SWH and "wsp" for wind speed
data_type = 'SIR_GOP_L2';
date0 = datenum('20140410','yyyymmdd');
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
            date0 = datenum(2014,4,10); % no CryoSat2 data prior to that date
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
            path1 = [path2Code num2str(Y) '_wsp_swh/'];
        else
            path1 = [path2Code num2str(Y) '_wsp_swh/',datestr(date0,'mmm'),'/'];
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
        wsp = cell(length(fn),1);
        swh = cell(length(fn),1);
        t = cell(length(fn),1);
        format1 = repmat('%d',[1 5]);
        format2 = repmat('%f',[1 13]);
        for i=1:length(fn)
            try
                fid = fopen([path1 fn{i}],'r'); %open data file
                header = fgetl(fid);
                header = regexp(header,'\s+','split');
                dum = fgetl(fid); %#ok discard 2nd row (contains units)
                [~,k,~] = intersect(header,{'WSPD','WVHT'}); %identify column
                data = fscanf(fid,[format1 format2]); %read data
                data = reshape(data,18,length(data)/18);
                tmp1 = data(k(1),:);
                tmp1(tmp1 > 95) = NaN;
                tmp2 = data(k(2),:);
                tmp2(tmp2 > 95) = NaN;
                wsp{i} = tmp1;
                swh{i} = tmp2;
                t{i} = ([data(1:5,:); zeros(1,size(data,2))])'; %add zeros for sec
                fclose(fid); %close file
            catch %#ok
                warning(['could not read data for buoy ',fn{i}]) %#ok
            end
        end
        
        kwsp = cellfun(@(x) any(~isnan(x)),wsp);
        kswh = cellfun(@(x) any(~isnan(x)),swh);
        
        wsp = wsp(kwsp);
        swh = swh(kswh);
        twsp = t(kwsp);
        tswh = t(kswh);
        lonwsp = lon(kwsp);
        latwsp = lat(kwsp);
        lonswh = lon(kswh);
        latswh = lat(kswh);
        
        % convert time to Matlab datenum
        twsp = cellfun(@datenum,twsp,'unif',0);
        tswh = cellfun(@datenum,tswh,'unif',0);
        
        
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
        
        % reject stations that are closer than "minDist" from coast
        dswh = distance2Coast(lonswh(:),latswh(:),xh(:),yh(:),1000)./1000;
        dwsp = distance2Coast(lonwsp(:),latwsp(:),xh(:),yh(:),1000)./1000;
        
        kswh = find(dswh > minDist);
        kwsp = find(dwsp > minDist);
        wsp = wsp(kwsp);
        swh = swh(kswh);
        twsp = twsp(kwsp);
        tswh = tswh(kswh);
        lonwsp = lonwsp(kwsp);
        latwsp = latwsp(kwsp);
        lonswh = lonswh(kswh);
        latswh = latswh(kswh);
        
        if strcmp(vari,'swh')
            V = swh;
            T = tswh;
            X = lonswh;
            Y = latswh;
        else
            V = wsp;
            T = twsp;
            X = lonwsp;
            Y = latwsp;
        end
        
        nstation = length(T);
        
        
        %         path_C2Data = '/scratch/general/cryosat/validation_data/';
        %         path3 = '/scratch/general/cryosat/daily_stats/'; %CryoSat-2 data
        
        
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
            ztg = V{ii}; %buoy data
            ttg = T{ii};
            xtg = X(ii);
            ytg = Y(ii);
            
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
                    xtemp = dataVal.(['x_' vari]);
                    ytemp = dataVal.(['y_' vari]);
                    vtemp = dataVal.(vari);
                    ttemp = dataVal.(['t_' vari]);
                    otemp = dataVal.(['orbit_' vari]);
                    
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
                    
                    kt = find(dtime < 1/48); % 30 minutes time difference
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
                o(k) = []; %#ok
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
        t_c2(kempty) = []; %#ok
        
        c2BuoyData.lon = lon;
        c2BuoyData.lat = lat;
        c2BuoyData.month = days{mm};
        c2BuoyData.v_c2 = v_c2;
        c2BuoyData.v_buoy = v_buoy; %#ok
        
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

K = find(v2 ~= 0 & v1 ~= 0 & ~(v1>2 & v2<0.2) & ~(v2>0.5 & v1>0.5 & v1./v2>3) &...
    v1<9 & v2<9);
v1 = v1(K);
v2 = v2(K);

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
lon = lon(K);
lat = lat(K);

% reject buoys located in the Great Lakes
K = find(~(lon>-93 & lon<-75 & lat>40 & lat<50));
lon = lon(K);
lat = lat(K);
v1 = v1(K);
v2 = v2(K);

if strcmp(vari,'swh')
    laby =  'CryoSat-2 SWH (m)';
    labx = 'Buoy SWH (m)';
    XYLIM = [0 10];
    xyTick = 0:2:20;
    units = ' m';
    xtxt = 5.5;
    ytxt = 1:0.75:3.25;
else
    laby = 'CryoSat-2 wind speed (m/s)';
    labx = 'Buoy wind speed (m)';
    XYLIM = [0 25];
    xyTick = 0:5:30;
    units = ' m/s';
    xtxt = 12;
    ytxt = 1.5:1.5:6;
end


% calculate OLS
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
    
    k = vtmp1~=0 & vtmp2 ~=0 & ~(lon>-93 & lon<-75 & lat>40 & lat<50) & ...
        ~(vtmp1>2 & vtmp2<0.2) & ~(vtmp2>0.5 & vtmp1>0.5 & vtmp1./vtmp2>3) & ...
        vtmp1<9 & vtmp2<9;
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
RMS = RMS.*100;

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
ylabel('RMS SWH (cm)')
set(gca,'FontSize',12,'XTick',xti,'XTickLabel',Month2,'YTick',0:dy:100)
ylim([ymin ymax])
rotateXLabels(gca(),45)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_' vari '_RMS'],'-r300') % save figure
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
ylabel('Slope SWH')
set(gca,'FontSize',12,'XTick',xti,'XTickLabel',Month2,'YTick',0:dy:2)
ylim([ymin ymax])
rotateXLabels(gca(),45)
%ylim([0.89 1.1])
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_' vari '_TAU'],'-r300') % save figure
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
ylabel(gca,laby,'FontSize',10)
xlabel(gca,labx,'FontSize',10)
set(gca,'LineWidth',1,'FontSize',10,'XTick',xyTick,'YTick',xyTick)
xlim(XYLIM)
ylim(XYLIM)
hold on
plot(tt,rols,'-','color',rgb2,'LineWidth',1.5)
plot(-1e2:1e2,-1e2:1e2,'--','color','black','LineWidth',1)
text(xtxt,ytxt(4),['R{^2} = ' R2])
text(xtxt,ytxt(3),['rms = ' rmsols units])
text(xtxt,ytxt(2),['slope = ' sols ' \pm ' esols])
text(xtxt,ytxt(1),['intcp = ' iols ' \pm ' eiols units])
box on
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_' vari '_ols'],'-r300') % save figure
close all

%plot for robustfit approach
rmsrob = sprintf('%4.2f',rmsrob);
srob = sprintf('%4.2f',brob(2));
irob = sprintf('%4.2f',brob(1));
esrob = sprintf('%4.2f',erob(2));
eirob = sprintf('%4.2f',erob(1));

figure(1)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
ha1 = axes('Units','centimeters','Position',[1.2 1 pos]);
scatter(v2,v1,15,rgb1,'filled');
ylabel(gca,laby,'FontSize',10)
xlabel(gca,labx,'FontSize',10)
set(gca,'LineWidth',1,'FontSize',10,'XTick',xyTick,'YTick',xyTick)
xlim(XYLIM)
ylim(XYLIM)
hold on
plot(tt,rrob,'-','color',rgb2,'LineWidth',1.5)
plot(-1e2:1e2,-1e2:1e2,'--','color','black','LineWidth',1)
text(xtxt,ytxt(4),['rms = ' rmsrob units])
text(xtxt,ytxt(3),['slope = ' srob ' \pm ' esrob])
text(xtxt,ytxt(2),['intcp = ' irob ' \pm ' eirob units])
box on

% inset plot
pa1 = get(ha1,'Position');
ha2 = axes('Units','centimeters','position', ...
    [pa1(1)+pa1(3)/10 pa1(2)+pa1(4)/1.65 pa1(3)/2.85 pa1(4)/2.85]);
scatter(ha2,v2,v1,15,rgb1,'filled');
set(ha2,'LineWidth',1,'FontSize',10,'XTick',0:0.5:1.5,'YTick',0:0.5:1.5,'box','on')
xlim(ha2,[0 1.5])
ylim(ha2,[0 1.5])
hold on
plot(ha2,tt,rrob,'-','color',rgb2,'LineWidth',1.5)
plot(ha2,-1e2:1e2,-1e2:1e2,'--','color','black','LineWidth',1)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_' vari '_robust'],'-r300') % save figure
close all

% plot differences
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
print('-depsc',[dir_out 'Fig_' vari '_diff'],'-r300') % save figure
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
print('-depsc',[dir_out 'Fig_buoy_' vari '_location'],'-r300') % save figure
close all
% -------------------------------------------------------------------------
