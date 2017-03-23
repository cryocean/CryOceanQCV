not sure this is still in use
function monthly_report_v1_nodisplay(data_type,report_date) 
addpath(genpath('/noc/users/fmc1q07/QCV_Cryo2/code/'))
%report_date = '20140913'; %yyyymmdd
%data_type = 'SIR_FDM_L2';

in_dir = ['/scratch/general/cryosat/' data_type '/'];
path1 = '/scratch/general/cryosat/daily_stats/';
path2 = '/scratch/general/cryosat/daily_data/';
dir_out = '/noc/users/fmc1q07/QCV_Cryo2/code/monthly_report/figures/';
%------------------ input data for statistics -----------------------------
path4 = '/noc/users/fmc1q07/QCV_Cryo2/code/';
load([path4 'DTU10MSS.mat']) % dtu10 MSS for the 20-Hz data
load([path4 'gshhs_i.mat']) % shoreline data
load([path4 'Groundtrack_20140101T000001_20150629T050000_2sec.mat']);
% -------------------------------------------------------------------------
prod = data_type(5:7);
fid = fopen([dir_out prod '_reportData.txt'],'w');
fid2 = fopen([dir_out prod '_warnings.txt'],'w');

%% Ensure data are available in matlab structure, if not process them
datef = datenum(report_date,'yyyymmdd');
date0 = datef -29; %30-day window
t_window = datestr(date0:datef,'yyyymmdd');
files = [repmat([data_type '_'],30,1) t_window repmat('.mat',30,1)];
kempty = false(1,length(files));
for i=1:length(files)
    if ~exist([path2 files(i,:)],'file')
        tic
        read_cryo2Data_v2(in_dir,path2,t_window(i,:),t_window(i,:));
        toc
    end
    if i == length(files) && ~exist([path2 files(i,:)],'file')
        error('dailyReport:noData','There are no data for present day. Exiting...') % exit if no data for present day
    elseif ~exist([path2 files(i,:)],'file')
        dailyEmptyStats(t_window(i,:),data_type,path1); % create empty structure if no data
        kempty(i) = true;
    end
end

%% Read statistics in 30-day window 
files = [repmat(['Stats_' data_type '_'],30,1) t_window repmat('.mat',30,1)];
data = struct([]);
for i=1:length(files) 
    if ~exist([path1 files(i,:)],'file')
        daily_stats_v5(t_window(i,:),data_type,path2,in_dir,path1, ...
                       DTU10MSS,gshhs_i,A);
    end
    data = [data load([path1 files(i,:)])]; %#ok
end

polar = data(end).polarMask;
orbit_day = data(end).orbits;
orbB = data(end).orbitBias_ssha;
orbitUni = unique(orbit_day);
orbf = data(end).orbitLastFull;
orb0 = data(end).orbitFirstFull;
orb1Day = orbit_day(1);
if orbf == 0
    orbit_day(orbit_day == orbit_day(end)) = [];
end
if orb0 == 0
    orbit_day(orbit_day == orbit_day(1)) = [];
end
perc_orbOcean = vertcat(data.perc_orbOcean);
perc_LRMOrb = vertcat(data.perc_LRMOrb);
perc_orbKeys = cell2mat(keys(perc_orbOcean));
perc_orbOcean = cell2mat(values(perc_orbOcean));
perc_LRMOrb = cell2mat(values(perc_LRMOrb));
data = rmfield(data,'perc_orbOcean');
data = rmfield(data,'perc_LRMOrb');
data = rmfield(data,'polarMask');
data = rmfield(data,'modeSeas');
flds = fieldnames(data);
data = cellfun(@(x) {horzcat(data.(x))},flds);
data = cell2struct(data,flds);
perc_LRMDay = data.perc_LRMDay;
recDate = date0:datef;
recDate(kempty) = [];
ndays = length(recDate);
ind = zeros(1,sum(data.nrec));
ind(cumsum([1 data.nrec(1:end-1)])) = 1;
ind = cumsum(ind);
recDate = recDate(ind);

% -------------------- prepare data for 2D histogram ----------------------
y_swh = data.swh;
y_noise2DHist_inLRM = data.noise_ssha4scatter_inLRM;
nan_2DHist_inLRM = data.nan2DHist_inLRM;
if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
    y_noise2DHist_outLRM = data.noise_ssha4scatter_outLRM;
    nan_2DHist_outLRM = data.nan2DHist_outLRM;
end
% -------------------------------------------------------------------------


%% Data latency
laten = data.laten;
tlaten = laten(1,:);
laten = laten(2,:);
if strcmp(data_type,'SIR_FDM_L2')
    units_laten = 'hours ';
    ytick = 0:0.5:30;
    xhtick = 0:1:9;
    xvalues = 0:0.25:9;
    xhist = [0 8];
    laten_thr = 3;
elseif strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2')
    laten = laten./24;
    units_laten = 'days ';
    if strcmp(data_type,'SIR_IOP_L2')
        ytick = 0.5:0.5:30;
        xvalues = 0:0.2:5;
        laten_thr = 2;
        xhist = [0.5 3.5];
        xhtick = 0:0.5:8;
    else
        ytick = 25:1:50;
        xvalues = 25:0.5:40;
        laten_thr = 30;
        xhist = [25 40];
        xhtick = 25:1:45;
    end
end

dates = cellstr(datestr(date0:datef,'dd/mm/yy')); % dates for the plot
%rgbn = linspecer(9,'qualitative');
rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
rgb3 = [220/255 221/255 216/255];

%plot box and whiskers
close all
tlatenU = unique(tlaten);
kj = tlaten == max(tlatenU);
mlaten = median(laten(kj)); % median of latency
latenup = mlaten+std(laten(kj)); % upper range median 1-sigma
latendown = mlaten-std(laten(kj)); % lower range median 1-sigma
fprintf(fid,'%3.1f %3.1f %3.1f\n',[mlaten latendown latenup]);
pr = zeros(2,length(tlatenU));
perc_3d = zeros(length(tlatenU),1); % percentage of data delivered in time
muLaten = zeros(length(tlatenU),1);
sigLaten = zeros(length(tlatenU),1);
ptile = [5 95]; % whiskers percentiles
for f=1:length(tlatenU)
    kj = tlaten == tlatenU(f);
    pr(:,f) = prctile(laten(kj),ptile);
    perc_3d(f) = sum(laten(kj) <= 3)/sum(kj)*100;
    muLaten(f) = mean(laten(kj));
    sigLaten(f) = std(laten(kj));
end
if perc_3d(end) < 100
    fprintf(fid2,'%s\n','latency_fail');
    fprintf(fid2,'%3.1f\n',perc_3d(end));
end
if muLaten(end) > laten_thr
    fprintf(fid2,'%s\n','latency_mean_high');
    fprintf(fid2,'%2.1f\n',muLaten(end));
end

pmax = max(pr(2,:));
pmin = min(pr(1,:));
xp = 1:length(muLaten);
if strcmp(data_type,'SIR_IOP_L2')
    yup = 3.7;
    yup2 = 8.3;
    ylab = '% delivered within 3 days';
elseif strcmp(data_type,'SIR_FDM_L2')
    yup = 5.2;
    yup2 = 12.3;
    ylab = '% delivered within 3 hours';
elseif strcmp(data_type,'SIR_GOP_L2')
    yup = 37;
    yup2 = 40;
    ylab = '% delivered within 3 days';
end


hold on
% --------------------------- compute ylim --------------------------------
yup3 = max(pmax+0.1*abs(pmax-pmin),yup);
ylim1 = max(pmin-0.1*abs(pmax-pmin),0.1);
ylim2 = min(yup2,yup3);
% -------------------------------------------------------------------------

% ---------------------- highlight present day box ------------------------
% Using patch resulted in strange behaviour after having added the
% threshold of 3 hours/days. As an alternative we use rectangle:
prect = [ndays-0.5,ylim1,1.5,ylim2-ylim1];
rectangle('Position',prect,'FaceColor',[1 1 0], ...
    'EdgeColor','none')
% -------------------------------------------------------- end highlighting
boxWhiskers(laten,tlaten,units_laten,dates,rgb3,'black',ptile)
set(findobj(gcf,'tag','Outliers'),'Visible','off') % remove outliers
ylim([ylim1 ylim2])
hs = shadedErrorBar(xp,muLaten',sigLaten',{'-','color',rgb2, ...
    'LineWidth',1.5,'markerfacecolor',rgb2},1);
set(gca,'YTick',ytick,'YAxisLocation','left','box','off')
set(hs.edge,'LineStyle','--','color',rgb2)
set(hs.patch,'Facecolor','none')
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 12.6 9.5]);
set(gca,'Units','centimeters','Position',[1.1 1.45 10.5 6.5])
hax1 = gca;
hax1_pos = get(hax1,'Position');
hax2 = axes('YAxisLocation','right','Color','none','LineWidth',1,...
    'YColor',rgb1,'XAxisLocation','top','XTick',0:5:ndays, ...
    'XTickLabel',[],'NextPlot','add','Xlim',[0 ndays+1],'Units','centimeters');
plot(hax2,xp,perc_3d,'Color',rgb1,'LineWidth',0.8)%right axis
set(hax2,'Position',hax1_pos,'Ylim',[min(perc_3d)-10 110],'YTick',0:5:110)
set(hax1,'YTick',0:0.5:10)
ytickl = get(hax1,'YTickLabel');
if size(ytickl,1) > 10
    set(hax1,'YTick',1:2:13);
end
% collocate the 100% tick in the right axis with the last tick in the
% left axis
l2 = min(perc_3d)-10;
lim1 = get(hax1,'ylim');
w1 = diff(lim1);
ytickl = str2num(get(hax1,'YTickLabel')); %#ok
if strcmp(data_type,'SIR_FDM_L2')
    d1 = lim1(2)-ytickl(end);
else
    d1 = lim1(2)-ytickl(end);
end
u2 = (100*w1-d1*l2)/(w1-d1);
set(hax2,'ylim',[l2 u2])
% ---------------------------- end collocating
ytickr = get(hax2,'YTickLabel');
if size(ytickr,1) > 10
    set(hax2,'YTick',0:20:110);
end
hlab = ylabel(ylab);
hlab_pos = get(hlab,'Position');
hlab_pos(1) = hlab_pos(1);
set(hlab,'Position',hlab_pos)
hlab2 = get(hax1,'ylabel');
hlab2_pos = get(hlab2,'Position');
hlab2_pos(1) = hlab2_pos(1)+0.2;
set(hlab2,'Position',hlab2_pos)
rotateXLabels(hax1,45)
% -----------------------------------highlight present day box with a patch
% yl = get(hax1,'ylim');
% xl = get(hax1,'xlim');
% xv = [29.5 29.5 xl(2) xl(2)];
% yv = [yl(1) yl(2) yl(2) yl(1)];
% patch(xv,yv,-1e+1.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
%     'EdgeColor','none','Parent',hax1);
% -------------------------------------------------------- end highlighting
plot(hax2,xp(end),perc_3d(end),'o','Color','black','MarkerSize',5, ...
    'MarkerFaceColor','black')%right axis
plot(hax1,[-1 50],[3 3],'Color','black','LineWidth',0.8)%left axis
set(hax1,'Layer','top')
set(gcf, 'PaperPositionMode', 'manual');
print('-depsc',[dir_out prod '_Fig_1'],'-r300') % save figure 1
close all

% ----------------- plot latency histogram for report day -----------------
figure(2)
kj = tlaten == max(tlatenU);
xx = laten(kj);
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 7.4 6.3]);
set(gca,'Units','centimeters','Position',[1.1 1.1 6 5])
xx(abs(xx) > 9) = NaN;
hist(xx,xvalues);
hh = findobj(gca,'Type','patch');
set(hh,'FaceColor',rgb2,'EdgeColor','w')
set(gca,'FontSize',8,'XTick',xhtick,'TickDir','out')
xlabel(['Latency (' units_laten ')'])
xlim(xhist)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out prod '_Fig_2'],'-r300') % save figure
close all
% ----------------------------------------------------end latency histogram


%% Data coverage and completeness
% total number of records for each day in last 30 days
nrec = data.nrec;

%clear t0f
fprintf(fid,'%-s\n',data.tStart{end});
fprintf(fid,'%-s\n',data.tStop{end});
fprintf(fid,'%u %u\n',[orbit_day(1) orbit_day(end)]);

% total number of records per orbit in last 30 days
orbitU = unique(data.orbits(1):data.orbits(end));
nrecOrbit = histc(data.orbits,orbitU);
nrecOrbitOcean = histc(data.orbits(data.surface_type == 0 | ...
    data.surface_type == 1),orbitU);
% remove first and last orbits if they are not complete
if data.orbitFirstFull(1) == 0
    orbitU(1) = [];
    nrecOrbit(1) = [];
    nrecOrbitOcean(1) = [];
end
if data.orbitLastFull(end) == 0
    orbitU(end) = [];
    nrecOrbit(end) = [];
    nrecOrbitOcean(end) = [];
end
[~,io] = intersect(perc_orbKeys,orbitU);
perc_orbOcean = perc_orbOcean(io);
perc_LRMOrb = perc_LRMOrb(io);
nth_orbit = 6324; %(5344 orbits/369 days and 1 record every 0.9434 s)
if strcmp(data_type,'SIR_FDM_L2')
    perc_orbit = nrecOrbit./(nth_orbit.*perc_LRMOrb./100).*100;
else
    perc_orbit = nrecOrbit./nth_orbit.*100; % percentage records per orbit
end
perc_orbitOcean = nrecOrbitOcean./(nth_orbit.*perc_orbOcean./100).*100;

korb = find(orbitU == orbitUni(1));
orb_thr = 80;

if any(perc_orbitOcean(korb:end) < orb_thr)
    porb = perc_orbitOcean(korb:end);
    orbU = orbitU(korb:end);
    k80 = find(porb < orb_thr);
    fprintf(fid2,'%s\n','orbit_dropout');
    for ko = 1:length(k80)
        fprintf(fid2,'%u %2.1f\n',[orbU(k80(ko)) porb(k80(ko))]);
    end
end
    
nth = 86400/0.9434; %theoretical number of records per day (1 record every 0.9434 s)
nth_ocean = nth.*data.perc_ocean./100; %theoretical over ocean

num_rec = nrec(end); % number of 1-Hz records for present day
if strcmp(data_type,'SIR_FDM_L2')
    perc_rec = num_rec/(nth*perc_LRMDay(end)/100)*100;
    nth = nth.*(perc_LRMDay./100);
else
    perc_rec = num_rec/nth*100; % precentage over theoretical max
end
nthOcean = round(nth_ocean(end));
nthTotal = round(nth(end));
fprintf(fid,'%u\n',nthTotal);
fprintf(fid,'%u\n',nthOcean);
fprintf(fid,'%u %4.1f\n',[num_rec perc_rec]);

nrec_ocean = data.nrec_ocean;
perc_recOcean = nrec_ocean./nth_ocean.*100; %percentage over theoretical max
fprintf(fid,'%u %4.1f\n',[nrec_ocean(end) perc_recOcean(end)]);

% generate warning if needed
ocean_thr = 80; %threshold
if perc_recOcean(end) < ocean_thr
    fprintf(fid2,'%s\n','ocean_dropout');
    fprintf(fid2,'%2.1f\n',perc_recOcean(end));
end


tt = (1:length(nrec))';
tr1 = regstats(nrec./nth.*100,tt);
e1 = sqrt(diag(tr1.covb));
e1 = e1(2);
tr1 = tr1.beta(2);
tr2 = regstats(perc_recOcean,tt);
e2 = sqrt(diag(tr2.covb));
e2 = e2(2);
tr2 = tr2.beta(2);
s1 = ['(',num2str(sprintf('%.2f',tr1)),' \pm ', ...
    num2str(sprintf('%.2f',e1)),' day^{-1})'];
s2 = ['(',num2str(sprintf('%.2f',tr2)),' \pm ', ...
    num2str(sprintf('%.2f',e2)),' day^{-1})'];
ylab = 'Percentage relative to theory (%) ';
leg = {['Total ' s1], ...
    ['Ocean/lake ' s2]};
pPos = [10.1 6.1];
pos = [8.7 5];

% plot percentage of records for each day in last 30-day window
plot_timeSeries(2,(1:ndays)',nrec./nth.*100,perc_recOcean,rgb1,rgb2,ylab, ...
    dates,leg,pPos,pos);
rotateXLabels(gca(),45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out prod '_Fig_3'],'-r300') % save figure 2
close all

% plot percentage of records per orbit in last 30-day window
[hp1,hp2] = plot_timeSeries(3,orbitU,perc_orbit,perc_orbitOcean,rgb1, ...
    rgb2,ylab,orbitU,leg,pPos,pos,'-');
delete([hp1 hp2])
%ylim([-15 100])
yl = get(gca,'ylim');
xl = get(gca,'xlim');
xv = [orb1Day orb1Day xl(2) xl(2)];
yv = [yl(1) yl(2) yl(2) yl(1)];
patch(xv,yv,-1e+9.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
    'EdgeColor','none','Parent',gca);
set(gca,'Layer','top')
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out prod '_Fig_4'],'-r300') % save figure 3
close all
% -------------------------------------- end data coverage and completeness


%% Diagnostics for SLA, SWH and SIGMA0
params = {data.ssha,data.swh,data.sigma0,data.wsp,data.misp};
%flag_block = data.flag_block;
vari = {'ssha','swh','sigma0','wsp','misp'};

% problems with orbit biases
if any(abs(orbB) > 100) % if mean of ssha over orbit is larger 0.7 m
    kbias = find(abs(orbB) > 100);
    fprintf(fid2,'%s\n','large_orbit_bias');
    for ko = 1:length(kbias)
        fprintf(fid2,'%u\n',orbitUni(kbias(ko)));
    end
end

cnt = 4; % counter for figures (three figure have already been saved)
for j=1:length(params) % loop over all parameters
    x = params{j};
    inLRM = data.inLRM;
    noPolar = data.noPolar;
    nan_x = data.(['validFlag_' vari{j}]);
    lonx = data.lon(nan_x);
    latx = data.lat(nan_x);
    recDatex = recDate(nan_x);
    if recDatex(end) ~= datef
        error('No record in the report date')
    end
    [~,ix] = unique(recDatex,'first');
    nrecx = data.(['nrec_' vari{j}]); % total number of valid 1-Hz
    if j == 5
        nrecPos = data.(['nrecPos_' vari{j}]); % total number of positive values
        percRecxTheorPos = (nrecPos./nth_ocean).*100;
    end
    num_recx = nrecx(end); % number of valid 1-Hz records for present day
    percRecx = (nrecx./nrec_ocean).*100; % percentage over #records ocean
    percRecxTheor = (nrecx./nth_ocean).*100; % percentage over theoret ocean
    
    % save data
    fprintf(fid,'%u\n',num_recx);
    fprintf(fid,'%4.1f\n',percRecx(end));
    fprintf(fid,'%4.1f\n',percRecxTheor(end));
    if j == 5
        fprintf(fid,'%4.1f\n',percRecxTheorPos(end));
    end
    
    % plot geographical distribution of valid data for present day
    cnt = cnt+1;
    units = [{'SSH anomaly (cm) '};{'SWH (m) '};{'Sigma0 (dB) '};...
        {'Altimeter wind speed (m/s)'};{'Square of off-nadir angle (10^{-2} deg^2) '}];
    titles = {'SSH anomaly','SWH','Sigma0','Altimeter wind speed','Square of off-nadir angle'};
    cytick = [5,1,1,2,1];
    scale = [100,1,1,1,100];
    xposTab = [4.5,5.7,5.4,5.5,5.5];
    
    strCol = {'p5','p25','median','p75','p95','mean','std'};
    nan_x2 = data.(['validNOCPolarIncluded_' vari{j}]);
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        % compute monthly statistics
        x1 = x(nan_x & nan_x2 & inLRM & noPolar);
        x2 = x(nan_x & nan_x2 & ~inLRM & noPolar);
        perctile2 = prctile(x2,[5 25 75 95]);
        
        cmax = [20,8,15,20,5];
        cmin = [-20,0,6,0,1];
        strCol = [{'Mode'} strCol]; %#ok
        dat_tab2 = {perctile2(1),perctile2(2),median(x2),perctile2(3), ...
            perctile2(4),mean(x2),std(x2)};
        dat_tab2 = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab2,'unif',false);
    else
        cmax = [10,8,15,20,5];
        cmin = [-30,0,6,0,1];
        x1 = x(nan_x & noPolar);
    end
    perctile1 = prctile(x1,[5 25 75 95]);
    dat_tab = {perctile1(1),perctile1(2),median(x1),perctile1(3), ...
            perctile1(4),mean(x1),std(x1)};
    
    dat_tab = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab,'unif',false);
    x = x(nan_x & nan_x2);
    yin = data.lat(nan_x & nan_x2);
    xin = data.lon(nan_x & nan_x2);
    zin = x.*scale(j); %
    if j == 5 && (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        xout = xin(zin < 0);
        yout = yin(zin < 0);
        xin = xin(zin > 0);
        yin = yin(zin > 0);
        zin = zin(zin > 0);
    end
    %titl = ['Flag-valid ' titles{j} ' ' data.date{end}];
    titl = [titles{j} ' measured by CryoSat-2, September 2014'];
    pPos = [13.9 9.2];
    pos = [10.5 9];
    %pPos = [15.3 10.1];
    %pos = [11.6 9.9];
    plot_map_nodisplay_monthly(cnt,xin,yin,zin,cmax(j),cmin(j),cytick(j), ...
        units{j},titl,pPos,pos);
    if j == 5 && (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        plotm(yout,xout,'o','MarkerSize',2,'Color',[0.3 0.3 0.3])
    end
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        xposTab = [4.8,5.7,5.4,5.5,5.5];
        sr = {'LRM','PLRM'};
        ht = plotTable(sr,strCol,[dat_tab;dat_tab2],[xposTab(j) 0.48],7,1);
        set(ht(:,1),'FontWeight','bold','BackGroundColor', ...
        [206/255 235/255 251/255])
        %plotm(polar(2,:),polar(1,:),'color',[0 0 0],'LineWidth',1)
    else
        ht = plotTable([],strCol,dat_tab,[xposTab(j) 0.48],8,1);
    end
    set(ht(1,:),'FontWeight','bold','BackGroundColor', ...
        [206/255 235/255 251/255])
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-painters','-r300') % save figure
    close all
    disp('paused')
    pause
    
    % plot Percentage over theoretical max records for each day
    cnt = cnt+1;
    units = [{'SSH anomaly '};{'SWH '};{'Sigma0 '};{'Wind speed'};...
        {'off-nadir angle '}]; %#ok
    
    tt = (1:length(percRecxTheor))';
    tr2 = regstats(percRecxTheor,tt);
    e2 = sqrt(diag(tr2.covb));
    e2 = e2(2);
    tr2 = tr2.beta(2);
    if j == 5
        tr1 = regstats(percRecxTheorPos,tt);
        e1 = sqrt(diag(tr1.covb));
        e1 = e1(2);
        tr1 = tr1.beta(2);
    end
    s1 = ['(',num2str(sprintf('%.2f',tr1)),' \pm ', ...
        num2str(sprintf('%.2f',e1)),' day^{-1})'];
    s2 = ['(',num2str(sprintf('%.2f',tr2)),' \pm ', ...
        num2str(sprintf('%.2f',e2)),' day^{-1})'];
    ylab = ['Percentage relative to theory (%) ']; %#ok
    if j == 5 && (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        leg = {['All values ' s2],['Only positive values ' s1]};
    else
        leg = ['Trend = ' s2(2:end-1)];
    end
    pPos = [10 6.1];
    pos = [8.7 5];
    if j == 5 && (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        plot_timeSeries(cnt,(1:ndays)',percRecxTheor,percRecxTheorPos,rgb2, ...
            rgb1,ylab,dates,leg,pPos,pos);
    else
        plot_percData(cnt,(1:ndays)',percRecxTheor,rgb2,ylab,dates,leg,pPos,pos);
    end
    rotateXLabels(gca(),45)
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
    close all
    
    % plot histogram
    cnt = cnt+1;
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        xtmp = params{j};
        nan_x1 = data.(['validFlag_' vari{j}]) & data.inLRM & data.noPolar;
        nan_x2 = data.(['validFlag_' vari{j}]) & ~data.inLRM & data.noPolar;
        x1 = xtmp(nan_x1);
        x2 = xtmp(nan_x2);
        recDatex1 = recDate(nan_x1);
        recDatex2 = recDate(nan_x2);
        [~,ix1] = unique(recDatex1,'first');
        [~,ix2] = unique(recDatex2,'first');
    else
        ix1 = ix;
        x1 = x;
        x2 = x;
        ix2 = ix;
    end
    
    
    xx = x1(ix1(end):end).*scale(j);
    xx2 = x2(ix2(end):end).*scale(j);
    if j == 5
        xx = xx(xx >= 0);
        xx2 = xx2(xx2 >= 0);
    end
    xlab = [{'SSH anomaly (cm)'};{'SWH (m)'};{'Sigma0 (dB)'}; ...
        {'Altimeter wind speed (m/s)'};{'Square of off-nadir angle (10^{-2} deg^2)'}];
    xmax = [50,12,20,24,6];
    xvalues = {-50:2:50,0:0.25:12,4:0.5:20,0:0.5:24,0:0.2:6};
    xmin = [-50,0,4,0,0];
    xticks = {xmin(1):10:xmax(1),xmin(2):2:xmax(2),xmin(3):2:xmax(3), ...
        xmin(4):2:xmax(4),xmin(5):1:xmax(5)};
    
    figure(cnt)
    set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 7.4 6.3]);
    set(gca,'Units','centimeters','Position',[1.1 1.1 6 5])
    switch j
        case 1
        xx(abs(xx) > 50) = NaN;
        xx2(abs(xx2) > 50) = NaN;
        case 2
        xx(xx > 12) = NaN;
        xx2(xx2 > 12) = NaN;
        case 3
        xx(xx > 20 | xx < 4) = NaN;
        xx2(xx2 > 20 | xx2 < 4) = NaN;
        case 4
        xx(xx > 24) = NaN;
        xx2(xx2 > 24) = NaN;
        case 5
        xx(xx > 6 | xx < xmin(j)) = NaN;
        xx2(xx2 > 6 | xx2 < xmin(j)) = NaN;
    end
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        [h1,x1] = hist(xx,xvalues{j});
        [h2,x2] = hist(xx2,xvalues{j});
        hold on
        hb1 = bar(x1,h1,'FaceColor',rgb2,'EdgeColor','none');
        hb2 = bar(x2,h2,0.5,'FaceColor',rgb1,'EdgeColor','none');
        box on
        hl = legend([hb1 hb2],'LRM','PLRM');
        set(hl,'box','off')
    else
        hist(xx,xvalues{j});
        hh = findobj(gca,'Type','patch');
        set(hh,'FaceColor',rgb2,'EdgeColor','w')
    end
    set(gca,'FontSize',8,'XTick',xticks{j},'TickDir','out')
    xlabel(xlab{j})
    xlim([xmin(j) xmax(j)])
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
    close all
    
    if j == 5 % mispointing
        cnt = cnt+1;
        units = [{'SSH anomaly (cm) '};{'SWH (m) '};{'Sigma0 (dB) '};...
            {'Altimeter wind speed (m/s)'};{'Square of off-nadir angle (10^{-2} deg^2) '}];
        ytick = {-50:10:50,0:2:20,0:2:20,0:2:20,-10:1:20};
        ylims = {[-50 40],[0 8],[7.5 16],[0 20],[0.5 5]};
        
        recDatex = recDatex(x >= 0);
        x = x(x >= 0);
        [recDatexU,~] = unique(recDatex,'first');
        pr = zeros(2,length(recDatexU));
        ptile = [5 95];
        for f=1:length(recDatexU)
            kj = recDatex == recDatexU(f);
            pr(:,f) = prctile(x(kj)*scale(j),ptile);
        end
        pmax = max(pr(2,:));
        pmin = min(pr(1,:));
        xp = 1:length(recDatexU);
        boxWhiskers(x*scale(j),recDatex,units{j},dates,rgb3,'black',ptile)
        hh = findobj(gcf,'tag','Outliers');
        set(hh,'Visible','off')
        hs = shadedErrorBar(xp,data.(['mean_' vari{j}]).*scale(j), ...
            data.(['std_' vari{j}]).*scale(j),{'-','color',rgb2,'LineWidth', ...
            1.5,'markerfacecolor',rgb2});
        set(hs.edge,'LineStyle','none')
        set(hs.patch,'Facecolor','none')
        ylim([pmin-0.1*abs(pmax-pmin) pmax+0.1*abs(pmax-pmin)])
        set(gca,'YTick',ytick{j})
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 10.3 6.7]);
        set(gca,'Units','centimeters','Position',[1.4 1.25 8.7 5])
        ylim(ylims{j})
        rotateXLabels(gca,45)
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        xv = [ndays-0.5 ndays-0.5 xl(2) xl(2)];
        yv = [yl(1) yl(2) yl(2) yl(1)];
        patch(xv,yv,-1e+9.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
            'EdgeColor','none','Parent',gca);
        set(gca,'Layer','top')
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
    end
    
    if(j < 5) % exclude mispointing
        nan_x = data.(['validFlag_' vari{j}]) & data.noPolar;
        lonx = data.lon(nan_x);
        latx = data.lat(nan_x);
        recDatex = recDate(nan_x);
        [~,ix] = unique(recDatex,'first');
        if j < 4
            scaleNoise = [100,100,100];
            noiseDiffNoc = 0;
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                recDatex1 = recDate(nan_x & data.inLRM);
                recDatex2 = recDate(nan_x & ~data.inLRM);
                [~,ix1] = unique(recDatex1,'first');
                [~,ix2] = unique(recDatex2,'first');
                stdX1 = (data.(['noise_inLRM_' vari{j}])).*scaleNoise(j);
                meanStdX1 = (data.(['noiseavg_inLRM_' vari{j}])).*scaleNoise(j);
                %noiseDiff1 = (data.(['noiseDiff_inLRM_' vari{j}])).*scaleNoise(j);
                stdX2 = (data.(['noise_outLRM_' vari{j}])).*scaleNoise(j);
                meanStdX2 = (data.(['noiseavg_outLRM_' vari{j}])).*scaleNoise(j);
                %noiseDiff2 = (data.(['noiseDiff_outLRM_' vari{j}])).*scaleNoise(j);
                
                meanStdXNoc1 = (data.(['noiseavgNOC_inLRM_' vari{j}])).*scaleNoise(j);
                noiseStdNoc1 = (data.(['noiseStdNOC_inLRM_' vari{j}])).*scaleNoise(j);
                meanStdXNoc2 = (data.(['noiseavgNOC_outLRM_' vari{j}])).*scaleNoise(j);
                noiseStdNoc2 = (data.(['noiseStdNOC_outLRM_' vari{j}])).*scaleNoise(j);
            end
            stdX = (data.(['noise_' vari{j}])).*scaleNoise(j);
            meanStdX = (data.(['noiseavg_' vari{j}])).*scaleNoise(j);
            %noiseStd = (data.(['noiseStd_' vari{j}])).*scaleNoise(j);
            noiseDiff = (data.(['noiseDiff_' vari{j}])).*scaleNoise(j);
            %noiseDiffNoc = (data.(['noiseDiffNOC_' vari{j}])).*scaleNoise(j);
            meanStdXNoc = (data.(['noiseavgNOC_' vari{j}])).*scaleNoise(j);
            noiseStdNoc = (data.(['noiseStdNOC_' vari{j}])).*scaleNoise(j);
            
            % save data
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                fprintf(fid,'%4.1f\n',meanStdX1(end));
                fprintf(fid,'%4.1f\n',meanStdX2(end));
                fprintf(fid,'%4.1f\n',meanStdX1(end)/sqrt(20));
                fprintf(fid,'%4.1f\n',meanStdX2(end)/sqrt(20));
                fprintf(fid,'%4.1f\n',noiseDiff(end));
                fprintf(fid,'%4.1f\n',meanStdXNoc1(end));
                fprintf(fid,'%4.1f\n',meanStdXNoc2(end));
                fprintf(fid,'%4.1f\n',meanStdXNoc1(end)/sqrt(20));
                fprintf(fid,'%4.1f\n',meanStdXNoc2(end)/sqrt(20));
                fprintf(fid,'%4.1f\n',noiseDiffNoc(end));
            else
                fprintf(fid,'%4.1f\n',meanStdX(end));
                fprintf(fid,'%4.1f\n',meanStdX(end)/sqrt(20));
                fprintf(fid,'%4.1f\n',noiseDiff(end));
                fprintf(fid,'%4.1f\n',meanStdXNoc(end));
                fprintf(fid,'%4.1f\n',meanStdXNoc(end)/sqrt(20));
                fprintf(fid,'%4.1f\n',noiseDiffNoc(end));
            end
            
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                meanStdXNoc = meanStdXNoc1;
                noiseStdNoc = noiseStdNoc1;
            end
        
            % plot geographical distribution of 20-Hz SLA noise for present day
            cnt = cnt+1;
            units = [{'SSH anomaly noise (cm)'};{'SWH noise (cm)'}; ...
                {'Sigma0 noise (10^{-2} dB)'}];
            titles = {'Flag-valid SSH anomaly noise','Flag-valid SWH noise', ...
                'Flag-valid sigma0 noise'};
            cytick = [2,10,2];
            zin = stdX(ix(end):end);
            zin1 = zin;
            xin = lonx(ix(end):end);
            yin = latx(ix(end):end);
            
            strCol = {'p5','p25','median','p75','p95','mean','std'};
            
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                cmax = [17,90,20];
                cmin = [4,30,6];
                zin1 = stdX1(ix1(end):end);
                zin2 = stdX2(ix2(end):end);
                strCol = [{'Mode'} strCol]; %#ok
                dat_tab2 = {prctile(zin2,5),prctile(zin2,25),prctile(zin2,50), ...
                    prctile(zin2,75),prctile(zin2,95),nanmean(zin2),nanstd(zin2)};
                dat_tab2 = cellfun(@(x) sprintf('%0.1f',x),dat_tab2,'unif',false);
            else
                cmax = [12,80,18];
                cmin = [4,30,6];
            end
            
            dat_tab = {prctile(zin1,5),prctile(zin1,25),prctile(zin1,50), ...
                prctile(zin1,75),prctile(zin1,95),nanmean(zin1),nanstd(zin1)};
            dat_tab = cellfun(@(x) sprintf('%0.1f',x),dat_tab,'unif',false);
            titl = [titles{j} ' ' data.date{end}];
            pPos = [13.9 9.2];
            pos = [10.5 9];
            plot_map_nodisplay(cnt,xin,yin,zin,cmax(j),cmin(j),cytick(j), ...
                units{j},titl,pPos,pos);
            
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                sr = {'LRM','PLRM'};
                ht = plotTable(sr,strCol,[dat_tab;dat_tab2],[5.0 0.48],7,1);
                set(ht(:,1),'FontWeight','bold','BackGroundColor', ...
                    [206/255 235/255 251/255])
                plotm(polar(2,:),polar(1,:),'color',[0 0 0],'LineWidth',1)
            else
                ht = plotTable([],strCol,dat_tab,[5.2 0.48],8,1);
            end
            
            set(ht(1,:),'FontWeight','bold','BackGroundColor', ...
                [206/255 235/255 251/255])
            set(gcf, 'PaperPositionMode', 'manual')
            print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
            close all
        end
        
        if j==1
            % plot 2D histogram of noise in SLA vs SWH
            y = y_swh; % values for present day only
            nan_xy = nan_2DHist_inLRM; % values for present day only
            y = y(nan_xy);
            ssh_noise = y_noise2DHist_inLRM; % values for present day only
            %ssh_noise(y == 0) = NaN;
            %y(y == 0) = NaN;
            cnt = cnt+1;
            figure(cnt);
            set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 11.5 9.5]);
            set(gca,'Units','centimeters','Position',[1.2 1.1 10 8])
            axis_hist = [0 10 0 20];
            z = [y' ssh_noise'.*100];
            [cx,xmedian] = hist2d(z,axis_hist,133,267);
            hold on
            plot(cx{1},xmedian,'LineWidth',2,'color','black')
            box on
            pu = get(gcf,'PaperUnits');
            pp = get(gcf,'PaperPosition');
            set(gcf,'Units',pu,'Position',pp)
            set(gca,'LineWidth',1.5,'FontSize',10)
            ylabel('SSH anomaly noise (cm)')
            xlabel('SWH (m)')
            %xlim([0 10]) %30-day window
            xlim([0 8]) % only 1 day
            %ylim([0 20]) % 30-day window
            ylim([1 20])
            set(gca,'XTick',0:2:16,'YTick',0:2:20)
            hc = colorbar;
            xpos = get(gca,'position');
            %set(hc,'YTick',0:100:10000) % 30-day window
            set(hc,'YTick',0:200:3000) % only 1 day
            set(hc,'FontSize',10)
            cpos = get(hc,'Position');
            cpos(3) = 0.6.*cpos(3);
            cpos(2) = cpos(2)+0.708;
            set(hc,'Position',cpos)
            set(gca,'position',xpos)
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                text(6.7,2,'LRM');
            end
            set(gcf, 'PaperPositionMode', 'manual')
            print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
            close all
            
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                % plot 2D histogram for Pseudo LRM
                y = y_swh; % values for present day only
                nan_xy = nan_2DHist_outLRM; % values for present day only
                y = y(nan_xy);
                ssh_noise = y_noise2DHist_outLRM; % values for present day only
                ssh_noise(y == 0) = NaN;
                y(y == 0) = NaN;
                figure(cnt);
                set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 11.5 9.5]);
                set(gca,'Units','centimeters','Position',[1.2 1.1 10 8])
                axis_hist = [0 10 0 20];
                z = [y' ssh_noise'.*100];
                [cx,xmedian] = hist2d(z,axis_hist,133,267);
                hold on
                plot(cx{1},xmedian,'LineWidth',2,'color','black')
                box on
                pu = get(gcf,'PaperUnits');
                pp = get(gcf,'PaperPosition');
                set(gcf,'Units',pu,'Position',pp)
                set(gca,'LineWidth',1.5,'FontSize',10)
                ylabel('SSH anomaly noise (cm)')
                xlabel('SWH (m)')
                %xlim([0 10]) %30-day window
                xlim([0 8]) % only 1 day
                %ylim([0 20]) % 30-day window
                ylim([1 20])
                set(gca,'XTick',0:2:16,'YTick',0:2:20)
                hc = colorbar;
                xpos = get(gca,'position');
                %set(hc,'YTick',0:100:10000) % 30-day window
                set(hc,'YTick',0:20:300) % only 1 day
                set(hc,'FontSize',10)
                cpos = get(hc,'Position');
                cpos(3) = 0.6.*cpos(3);
                cpos(2) = cpos(2)+0.708;
                set(hc,'Position',cpos)
                set(gca,'position',xpos)
                text(6.7,2,'PLRM')
                set(gcf, 'PaperPositionMode', 'manual')
                print('-depsc',[dir_out prod '_Fig_9b'],'-r300') % save figure
                close all
            end
            disp('stop')
            pause
        end
        
        % ---------------------- Validity by NOC --------------------------
        nan_x2 = data.(['validNOC_' vari{j}]);
        x = params{j};
        lonout = data.lon(~nan_x2 & nan_x);
        latout = data.lat(~nan_x2 & nan_x);
        recDateOut = recDate(~nan_x2 & nan_x);
        nan_x2 = nan_x2 & nan_x;
        x = x(nan_x2);
        lonx = data.lon(nan_x2);
        latx = data.lat(nan_x2);
        recDatex = recDate(nan_x2);
        [recDatexU,ix] = unique(recDatex,'first');
        [~,ixout] = unique(recDateOut,'first');
        nrecx = data.(['nrecNOC_' vari{j}]); % total number of valid 1-Hz
        num_recx = nrecx(end); % number of valid 1-Hz records for present day
        percRecx = (nrecx./nrec_ocean).*100; % percentage over #records ocean
        percRecxTheor = (nrecx./nth_ocean).*100; % percentage over theoret ocean
        % save data
        fprintf(fid,'%u\n',num_recx);
        fprintf(fid,'%4.1f\n',percRecx(end));
        fprintf(fid,'%4.1f\n',percRecxTheor(end));
        
        
        % plot geographical distribution of valid data for present day
        cnt = cnt+1;
        units = [{'SSH anomaly (cm) '};{'SWH (m) '};{'Sigma0 (dB) '};...
        {'Altimeter wind speed (m/s)'}];
        titles = {'SSH anomaly','SWH','Sigma0','Altimeter wind speed'};
        cytick = [5,1,1,2];
        scale = [100,1,1,1];
        xposTab = [4.5,5.7,5.4,5.5];
        
        strCol = {'p5','p25','median','p75','p95','mean','std'};
        
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            cmax = [20,8,14,20];
            cmin = [-20,0,7,0];
            strCol = [{'Mode'} strCol]; %#ok
            d1 = 'inLRM_';
            d2 = 'outLRM_';
            dat_tab2 = {data.(['p5NOC_' d2 vari{j}])(end),data.(['p25NOC_' d2 vari{j}])(end) ...
                ,data.(['medianNOC_' d2 vari{j}])(end),data.(['p75NOC_' d2 vari{j}])(end), ...
                data.(['p95NOC_' d2 vari{j}])(end),data.(['meanNOC_' d2 vari{j}])(end), ...
                data.(['stdNOC_' d2 vari{j}])(end)};
            dat_tab2 = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab2,'unif',false);
        else
            cmax = [10,8,14,20];
            cmin = [-30,0,7,0];
            d1 = '';
        end
        
        dat_tab = {data.(['p5NOC_' d1 vari{j}])(end),data.(['p25NOC_' d1 vari{j}])(end) ...
            ,data.(['medianNOC_' d1 vari{j}])(end),data.(['p75NOC_' d1 vari{j}])(end), ...
            data.(['p95NOC_' d1 vari{j}])(end),data.(['meanNOC_' d1 vari{j}])(end), ...
            data.(['stdNOC_' d1 vari{j}])(end)};
        dat_tab = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab,'unif',false);
        yin = latx(ix(end):end);
        xin = lonx(ix(end):end);
        zin = x(ix(end):end).*scale(j); %
        titl = ['Science-valid ' titles{j} ' ' data.date{end}];
        pPos = [13.9 9.2];
        pos = [10.5 9];
        plot_map_nodisplay(cnt,xin,yin,zin,cmax(j),cmin(j),cytick(j), ...
            units{j},titl,pPos,pos);
        if ~isempty(ixout)
            plotm(latout(ixout(end):end),lonout(ixout(end):end),'o', ...
                'MarkerSize',2.5,'Color',[0.4 0.4 0.4])
        end
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            xposTab = [4.8,5.7,5.4,5.5,5.5];
            sr = {'LRM','PLRM'};
            ht = plotTable(sr,strCol,[dat_tab;dat_tab2],[xposTab(j) 0.48],7,1);
            set(ht(:,1),'FontWeight','bold','BackGroundColor', ...
                [206/255 235/255 251/255])
            plotm(polar(2,:),polar(1,:),'color',[0 0 0],'LineWidth',1)
        else
            ht = plotTable([],strCol,dat_tab,[xposTab(j) 0.48],8,1);
        end
        set(ht(1,:),'FontWeight','bold','BackGroundColor', ...
            [206/255 235/255 251/255])
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
        
        % plot histogram
        cnt = cnt+1;
        
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            xtmp = params{j};
            nan_x1 = nan_x2 & data.inLRM;
            nan_X2 = nan_x2 & ~data.inLRM;
            x1 = xtmp(nan_x1);
            x2 = xtmp(nan_X2);
            recDatex1 = recDate(nan_x1);
            recDatex2 = recDate(nan_X2);
            [~,ix1] = unique(recDatex1,'first');
            [~,ix2] = unique(recDatex2,'first');
        else
            ix1 = ix;
            x1 = x;
            x2 = x;
            ix2 = ix;
        end
        
        xx = x1(ix1(end):end).*scale(j);
        xx2 = x2(ix2(end):end).*scale(j);
        
        xlab = [{'SSH anomaly (cm)'};{'SWH (m)'};{'Sigma0 (dB)'};{'Altimeter wind speed (m/s)'}];
        xvalues = {-50:2:50,0:0.25:12,4:0.5:20,0:0.5:24};
        xmin = [-50,0,4,0];
        xmax = [50,12,20,24];
        xticks = {xmin(1):10:xmax(1),xmin(2):2:xmax(2),xmin(3):2:xmax(3), ...
            xmin(4):2:xmax(4)};
        
        figure(cnt)
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 7.4 6.3]);
        set(gca,'Units','centimeters','Position',[1.1 1.1 6 5])
        switch j
            case 1
                xx(abs(xx) > 50) = NaN;
                xx2(abs(xx2) > 50) = NaN;
            case 2
                xx(xx > 12) = NaN;
                xx2(xx2 > 12) = NaN;
            case 3
                xx(xx > 20 | xx < 4) = NaN;
                xx2(xx2 > 20 | xx2 < 4) = NaN;
            case 4
                xx(xx > 24) = NaN;
                xx2(xx2 > 24) = NaN;
        end
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            [h1,x1] = hist(xx,xvalues{j});
            [h2,x2] = hist(xx2,xvalues{j});
            hold on
            hb1 = bar(x1,h1,'FaceColor',rgb2,'EdgeColor','none');
            hb2 = bar(x2,h2,0.5,'FaceColor',rgb1,'EdgeColor','none');
            box on
            hl = legend([hb1 hb2],'LRM','PLRM');
            set(hl,'box','off')
        else
            hist(xx,xvalues{j});
            hh = findobj(gca,'Type','patch');
            set(hh,'FaceColor',rgb2,'EdgeColor','w')
        end
        set(gca,'FontSize',8,'XTick',xticks{j},'TickDir','out')
        xlabel(xlab{j})
        xlim([xmin(j) xmax(j)])
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
        
        
        % plot box whiskers
        cnt = cnt+1;
        units = [{'SSH anomaly (cm) '};{'SWH (m) '};{'Sigma0 (dB) '};...
        {'Altimeter wind speed (m/s)'}];
        ytick = {-50:10:50,0:2:20,0:2:20,0:2:20};
        ylims = {[-30 30],[0 8],[7.5 16],[0 20]};
        pr = zeros(2,length(recDatexU));
        ptile = [5 95];
        for f=1:length(recDatexU)
            kj = recDatex == recDatexU(f);
            pr(:,f) = prctile(x(kj)*scale(j),ptile);
        end
        pmax = max(pr(2,:));
        pmin = min(pr(1,:));
        xp = 1:length(recDatexU);
        if j ~= 1
            boxWhiskers(x*scale(j),recDatex,units{j},dates,rgb3,'black',ptile)
            hh = findobj(gcf,'tag','Outliers');
            set(hh,'Visible','off')
            hs = shadedErrorBar(xp,data.(['meanNOC_' vari{j}]).*scale(j), ...
                data.(['stdNOC_' vari{j}]).*scale(j),{'-','color',rgb2,'LineWidth', ...
                1.5,'markerfacecolor',rgb2});
            set(hs.edge,'LineStyle','none')
            set(hs.patch,'Facecolor','none')
            
        else
            if strcmp(data_type,'SIR_FDM_L2')
                ylims = {[-40 20],[0 8],[7.5 16],[0 20]};
            end
            xv = [ndays-0.5   ndays-0.5   ndays+1   ndays+1];
            yv = [ylims{j}    fliplr(ylims{j})];
            hold on
            patch(xv,yv,-1e+9.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
                'EdgeColor','none','Parent',gca);
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                hs1 = errorbar(xp,data.(['meanNOC_inLRM_' vari{j}]).*scale(j), ...
                    data.(['stdNOC_inLRM_' vari{j}]).*scale(j),'o','color',rgb2, ...
                    'LineWidth',0.7,'MarkerSize',4,'MarkerFaceColor',rgb2);
            else
                errorbar(xp,data.(['meanNOC_' vari{j}]).*scale(j), ...
                    data.(['stdNOC_' vari{j}]).*scale(j),'o','color',rgb2, ...
                    'LineWidth',0.7,'MarkerSize',4,'MarkerFaceColor',rgb2);
            end
            ylabel(units{j},'FontSize',10)
            set(gca,'LineWidth',1,'FontSize',10,'XTick',[1:5:ndays ndays], ...
                'XTickLabel',[dates(1:5:end);dates(end)])
            xlim([0 ndays+1])
            box on
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                hold on
                hs2 = shadedErrorBar(xp,data.(['meanNOC_outLRM_' vari{j}]).*scale(j), ...
                    data.(['stdNOC_outLRM_' vari{j}]).*scale(j),{'-o','color',rgb1,'LineWidth', ...
                    1.5,'markerfacecolor',rgb1,'MarkerSize',3});
                set(hs2.edge,'LineStyle','-','color',rgb1,'LineWidth',0.5)
                set(hs2.patch,'Facecolor','none')
                hl = legend([hs1,hs2.mainLine],'LRM','PLRM','Location','northWest'); %#ok
                %set(hl,'box','off')
            end
        end
        ylim([pmin-0.1*abs(pmax-pmin) pmax+0.1*abs(pmax-pmin)])
        set(gca,'YTick',ytick{j})
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 10.3 6.7]);
        set(gca,'Units','centimeters','Position',[1.4 1.25 8.7 5])
        ylim(ylims{j})
        rotateXLabels(gca,45)
        if j ~= 1
            yl = get(gca,'ylim');
            xl = get(gca,'xlim');
            xv = [ndays-0.5 ndays-0.5 xl(2) xl(2)];
            yv = [yl(1) yl(2) yl(2) yl(1)];
            patch(xv,yv,-1e+9.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
                'EdgeColor','none','Parent',gca);
        end
        set(gca,'Layer','top')
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
        if j < 4
            % plot noise level for each day in last 30-day window
            cnt = cnt+1;
            units = [{'SSH anomaly noise (cm)'};{'SWH noise (cm)'}; ...
                {'Sigma0 noise (10^{-2} dB)'}];
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                ymax = max(max(meanStdXNoc),max(meanStdXNoc2)) ...
                    +max(max(noiseStdNoc),max(noiseStdNoc2));
                ymin = min(min(meanStdXNoc),min(meanStdXNoc2)) ...
                    -max(max(noiseStdNoc),max(noiseStdNoc2));
                dy = ymax-ymin;
                ylims = [ymin-0.1*dy ymax+0.4*dy];
            else
                ymax = max(meanStdXNoc)+max(noiseStdNoc);
                ymin = min(meanStdXNoc)-max(noiseStdNoc);
                dy = ymax-ymin;
                ylims = [ymin-0.1*dy ymax+0.1*dy];
            end
            ytick = {4:0.2:20,30:1:80,0:0.5:20};
            ytick2 = {4:2.0:20,25:10:160,0:2:26};
            xv = [ndays-0.5   ndays-0.5   ndays+1   ndays+1];
            yv = [ylims    fliplr(ylims)];
            hold on
            patch(xv,yv,-1e+9.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
                'EdgeColor','none','Parent',gca);
            hs1 = errorbar((1:ndays)',meanStdXNoc,noiseStdNoc,'o','color',rgb2, ...
                'LineWidth',0.7,'MarkerSize',4,'MarkerFaceColor',rgb2);
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                hold on
                hs2 = shadedErrorBar(xp,meanStdXNoc2,noiseStdNoc2, ...
                    {'-o','color',rgb1,'LineWidth', ...
                    1.5,'markerfacecolor',rgb1,'MarkerSize',3});
                set(hs2.edge,'LineStyle','-','color',rgb1,'LineWidth',0.5)
                set(hs2.patch,'Facecolor','none')
                hl = legend([hs1,hs2.mainLine],'LRM','PLRM','Location','northWest'); %#ok
            end
            ylabel(units{j},'FontSize',10)
            set(gca,'LineWidth',1,'FontSize',10,'XTick',[1:5:ndays ndays], ...
                'XTickLabel',[dates(1:5:end);dates(end)])
            xlim([0 ndays+1])
            box on
            set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 10.3 6.7]);
            set(gca,'Units','centimeters','Position',[1.4 1.25 8.7 5])
            ylim(ylims)
            set(gca,'YTick',ytick{j})
            ytickl = get(gca,'YTickLabel');
            if size(ytickl,1) > 10
                set(gca,'YTick',ytick2{j});
            end
            rotateXLabels(gca(),45)
            set(gca,'Layer','top')
            set(gcf, 'PaperPositionMode', 'manual')
            print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
            close all
        end
        
        % crossover analysis for SSH
        if j == 1
            %diffx = data.diffCross;
            %xcross = data.lonCross;
            %ycross = data.latCross;
            meanx = data.meanCross.*scaleNoise(j);
            stdx = data.stdCross.*scaleNoise(j);
            cnt = cnt+1;
            units = [{'Mean crossovers (cm)'};{'SWH noise (cm)'}; ...
                {'Sigma0 noise (10^{-2} dB)'}];
            ymax = max(meanx)+max(stdx);
            ymin = min(meanx)-max(stdx);
            dy = ymax-ymin;
            %ylims = [ymin(j)-0.1*dy(j) ymax(j)+0.1*dy(j)];
            ylims = [0 ymax+0.1*dy]; % absolute values are >0
            ytick = {0:2:30,30:1:80,0:0.3:20};
            ytick2 = {0:5:50,30:2:80,0:0.5:20};
            xv = [ndays-0.5   ndays-0.5   ndays+1   ndays+1];
            yv = [ylims    fliplr(ylims)];
            hold on
            patch(xv,yv,-1e+9.*ones(size(xv)),[1 1 0],'FaceAlpha',1, ...
                'EdgeColor','none','Parent',gca);
            errorbar((1:ndays)',meanx,stdx,'o','color',rgb2, ...
                'LineWidth',0.7,'MarkerSize',4,'MarkerFaceColor',rgb2);
            ylabel(units{j},'FontSize',10)
            set(gca,'LineWidth',1,'FontSize',10,'XTick',[1:5:30 30], ...
                'XTickLabel',[dates(1:5:end);dates(end)])
            xlim([0 ndays+1])
            box on
            set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 10.3 6.7]);
            set(gca,'Units','centimeters','Position',[1.4 1.25 8.7 5])
            ylim(ylims)
            set(gca,'YTick',ytick{j})
            ytickl = get(gca,'YTickLabel');
            if size(ytickl,1) > 10
                set(gca,'YTick',ytick2{j});
            end
            rotateXLabels(gca(),45)
            set(gca,'Layer','top')
            set(gcf, 'PaperPositionMode', 'manual')
            print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
            close all
        end
        % --------------------------------------------- end validity by NOC
    end
    
end
fclose(fid);
fclose(fid2);
% ------------------------------------------------- end diagnostics


