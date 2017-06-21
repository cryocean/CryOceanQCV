function monthly_report_IOP_FDM_v3_nodisplay_TDSv2(data_type,report_date)

% TDS working version as of 7/6/17 - CB
% ADDED THIS FOR TDS
addpath /noc/mpoc/smos/software/tools_cb/


addpath(genpath('/noc/users/cryo/QCV_Cryo2/code/'))
% report_date = 'Jul-2013'; %mmm-yyyy (month-year)
% data_type='SIR_GOP_L2';
% UPDATE HERE JAN 2017
% in_dir = ['/scratch/general/cryosat/' data_type '/'];
% path1 = '/scratch/general/cryosat/daily_stats/';
% path2 = '/scratch/general/cryosat/daily_data/';
% path3 = '/scratch/general/cryosat/crossOverMonthly/';

in_dir = ['/noc/mpoc/cryo/TDS_March_2017/' data_type '/'];
path1 = '/noc/mpoc/cryo/cryosat/daily_stats/TDS/';
warning('possibly change these two lines')
path2 = '/noc/mpoc/cryo/TDS_testout/';
path3 = '/noc/mpoc/cryo/cryosat/crossOverMonthly/';
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/TDS/';

date0 = datenum(report_date,'mmm-yyyy');
[Y,M] = datevec(date0);
datef = date0+eomday(Y,M)-1;
t_window = datestr(date0:datef,'yyyymmdd');
ndays = size(t_window,1); % number of days in month

%------------------ input data for statistics -----------------------------
path4 = '/noc/users/cryo/QCV_Cryo2/code/';
load([path4 'DTU10MSS.mat']) % dtu10 MSS for the 20-Hz data
load([path4 'gshhs_i.mat']) % shoreline data

t_window2 = datestr((date0-11):(datef+11),'yyyymmdd');
fileTracks = [repmat('groundTrack_',ndays+22,1) t_window2 repmat('.mat',ndays+22,1)];
A = [];
for i=1:size(fileTracks,1)
    a = load([path4 'groundTracks/' fileTracks(i,:)]);
    A = [A a.A]; %#ok
end
% -------------------------------------------------------------------------

prod = data_type(5:7);
fid = fopen([dir_out prod '_reportData.txt'],'w');
fid2 = fopen([dir_out prod '_warnings.txt'],'w');

%% delete unused matlab structure files
dateEnd = date0-30;
dateInit = datenum(2014,04,10);
t_windowDel = datestr(dateInit:dateEnd,'yyyymmdd');
ndel = size(t_windowDel,1);
files = [repmat([data_type '_'],ndel,1) t_windowDel repmat('.mat',ndel,1)];
% for i=1:size(files,1)
%     if exist([path2 files(i,:)],'file')
%         delete([path2 files(i,:)])
%         delete([path2 'sph_' files(i,:)])
%         delete([path2 'mph_' files(i,:)])
%     end
% end

%% Ensure data are available in matlab structure, if not process them
files = [repmat([data_type '_'],ndays,1) t_window repmat('.mat',ndays,1)];
files_stat = [repmat(['Stats_' data_type '_'],ndays,1) t_window repmat('.mat',ndays,1)];
kempty = false(1,length(files));
for i=1:size(files,1)
    if ~exist([path2 files(i,:)],'file') && ~exist([path1 files_stat(i,:)],'file')
        read_cryo2DatavC(in_dir,path2,t_window(i,:),t_window(i,:));
    end
    if ~exist([path2 files(i,:)],'file') && ~exist([path1 files_stat(i,:)],'file')
        dailyEmptyStats(t_window(i,:),data_type,path1); % create empty structure if no data
        kempty(i) = true;
    end
end
%
% got_to_this_point = 1 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point


%% Read statistics for month
data = struct([]);
corrections = struct([]);
for i=1:size(files_stat,1)
    %  load testdataweds.mat
    if ~exist([path1 files_stat(i,:)],'file')
        %         load fri_temp.mat
        daily_stats_v6_TDS(t_window(i,:),data_type,path2,in_dir,path1, ...
            DTU10MSS,gshhs_i,A);
    end
    data = [data load([path1 files_stat(i,:)])]; %#ok %this only works if there is no struct in struct
    tmp = load([path1 files_stat(i,:)]);
    corrections = [corrections tmp.corrections]; %#ok
    if strcmp(data_type,'SIR_GOP_L2')
        spectra.klrm{i} = tmp.inLRM;
        spectra.h_hi{i} = tmp.hi_ssha;
        spectra.t_hi{i} = tmp.hi_utc;
        spectra.f1{i} = tmp.validFlag_ssha;
        spectra.f2{i} = tmp.validNOC_ssha;
    end
end

polar = data(end).polarMask;
orbit_day = data(end).orbits;
orbB = data(end).orbitBias_ssha;
orbitUni = unique(orbit_day);
orbit0 = unique(data(1).orbits);
orbitf = unique(data(end).orbits);
orbf = data(end).orbitLastFull;
orb0 = data(1).orbitFirstFull;
if orbf == 0
    orbitf(end) = [];
end
if orb0 == 0
    orbit0(1) = [];
end


perc_orbOcean = vertcat(data.perc_orbOcean);
perc_LRMOrb = vertcat(data.perc_LRMOrb);
perc_outSARinOrb = vertcat(data.perc_outSARinOrb);
perc_orbKeys = cell2mat(keys(perc_orbOcean));
perc_orbOcean = cell2mat(values(perc_orbOcean));
perc_LRMOrb = cell2mat(values(perc_LRMOrb));
perc_outSARinOrb = cell2mat(values(perc_outSARinOrb));
data = rmfield(data,'perc_orbOcean');
data = rmfield(data,'perc_LRMOrb');
data = rmfield(data,'polarMask');
data = rmfield(data,'perc_outSARinOrb');
data = rmfield(data,'modeSeas');
data = rmfield(data,'corrections');
flds = fieldnames(data);
data = cellfun(@(x) {horzcat(data.(x))},flds);
data = cell2struct(data,flds);
flds = fieldnames(corrections);
corrections = cellfun(@(x) {horzcat(corrections.(x))},flds);
corrections = cell2struct(corrections,flds);
perc_LRMDay = data.perc_LRMDay;
perc_outSARin = data.perc_outSARin;
recDate = date0:datef;
recDate(kempty) = [];
ndays = length(recDate);


%% Data latency
laten = data.laten;
tlaten = laten(1,:);

laten = laten(1,:);
if strcmp(data_type,'SIR_FDM_L2')
    units_laten = 'hours';
    ytick = 0:0.5:30;
    xhtick = 0:1:14;
    xvalues = 0:0.25:14;
    xhist1 = [0 8];
    xhist2 = [0 14];
    laten_thr = 3;
    thr1 = 3;
elseif strcmp(data_type,'SIR_IOP_L2')
    laten = laten./24;
    units_laten = 'days';
    xvalues = 0:0.2:5;
    laten_thr = 2;
    xhist = [0.5 4];
    xhtick = 0:0.5:8;
    ytick = 0.5:0.5:30;
    thr1 = 3;
elseif strcmp(data_type,'SIR_GOP_L2')
    if  sum(~isnan(laten)) == 0 ; % ALL NAN
        all_latent_nan = true ;
    else
        laten = laten./24;
        units_laten = 'days';
        xvalues = 25:0.5:37;
        laten_thr = 30;
        xhist = [26 36];
        xhtick = 25:1:40;
        ytick = 24:1:50;
        thr1 = 30;
    end
end
% got_to_this_point = 2 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point

dates = cellstr(datestr(date0:datef,'dd/mm/yy')); % dates for the plot
%rgbn = linspecer(9,'qualitative');
rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
rgb3 = [220/255 221/255 216/255];

%plot box and whiskers
close all
if ~all_latent_nan
    tlatenU = unique(tlaten);
    mlaten = median(laten); % median of latency
    %latenup = mlaten+std(laten(kj)); % upper range median 1-sigma
    %latendown = mlaten-std(laten(kj)); % lower range median 1-sigma
    %latendown = max(latendown,0);
    latendown = min(laten);
    latenup = max(laten);
    fprintf(fid,'%3.1f %3.1f %3.1f\n',[mlaten latendown latenup]);
    pr = zeros(2,length(tlatenU));
    perc_3d = zeros(length(tlatenU),1); % percentage of data delivered in time
    muLaten = zeros(length(tlatenU),1);
    sigLaten = zeros(length(tlatenU),1);
    ptile = [5 95]; % whiskers percentiles
    for f=1:length(tlatenU)
        kj = tlaten == tlatenU(f);
        pr(:,f) = prctile(laten(kj),ptile); if isnan(pr(:,f));pr(:,f) = 0;end;
        perc_3d(f) = sum(laten(kj) <= thr1)/sum(kj)*100;if isnan(perc_3d(f));perc_3d(f) = 0;end;
        muLaten(f) = mean(laten(kj));if isnan(muLaten(f));muLaten(f) = 0;end;
        sigLaten(f) = std(laten(kj));if isnan(sigLaten(f));sigLaten(f) = 0;end;
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
        yup = 35.3;
        yup2 = 38.3;
        ylab = '% delivered within 30 days';
    end
    %     got_to_this_point = 2 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point
    
    
    % hold on
    % % --------------------------- compute ylim --------------------------------
    % yup3 = max(pmax+0.1*abs(pmax-pmin),yup);
    % if strcmp(data_type,'SIR_GOP_L2')
    %     ylim1 = max(pmin-0.1*abs(pmax-pmin),24);
    % else
    %     ylim1 = max(pmin-0.1*abs(pmax-pmin),0.1);
    % end
    % ylim2 = min(yup2,yup3);
    % % -------------------------------------------------------------------------
    %
    % boxWhiskers(laten,tlaten,units_laten,dates,rgb3,'black',ptile)
    % set(findobj(gcf,'tag','Outliers'),'Visible','off') % remove outliers
    % ylim([ylim1 ylim2])
    % hs = shadedErrorBar(xp,muLaten',sigLaten',{'-','color',rgb2, ...
    %     'LineWidth',1.5,'markerfacecolor',rgb2},1);
    % set(gca,'YTick',ytick,'YAxisLocation','left','box','off')
    % set(hs.edge,'LineStyle','--','color',rgb2)
    % set(hs.patch,'Facecolor','none')
    % set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 12.6 9.5]);
    % set(gca,'Units','centimeters','Position',[1 1.45 10.5 6.5])
    % hax1 = gca;
    % hax1_pos = get(hax1,'Position');
    % hax2 = axes('YAxisLocation','right','Color','none','LineWidth',1,...
    %     'YColor',rgb1,'XAxisLocation','top','XTick',0:5:ndays, ...
    %     'XTickLabel',[],'NextPlot','add','Xlim',[0 ndays+1],'Units','centimeters');
    % plot(hax2,xp,perc_3d,'Color',rgb1,'LineWidth',0.8)%right axis
    % set(hax2,'Position',hax1_pos,'Ylim',[min(perc_3d)-10 110],'YTick',0:5:110)
    % if strcmp(data_type,'SIR_GOP_L2')
    %     set(hax1,'YTick',24:0.5:40)
    % else
    %     set(hax1,'YTick',0:0.5:10)
    % end
    % ytickl = get(hax1,'YTickLabel');
    % if size(ytickl,1) > 10
    %     set(hax1,'YTick',1:2:13);
    %     if strcmp(data_type,'SIR_GOP_L2')
    %         set(hax1,'YTick',24:1:40);
    %     end
    % end
    %     got_to_this_point = 2.1 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point
    
    
    % collocate the 100% tick in the right axis with the last tick in the
    % left axis
    l2 = min(perc_3d)-10;
    lim1 = get(hax1,'ylim');
    w1 = diff(lim1);
    ytickl = str2num(get(hax1,'YTickLabel')); %#ok
    d1 = lim1(2)-ytickl(end);
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
    plot(hax1,[-1 50],[thr1 thr1],'Color','black','LineWidth',0.8)%left axis
    set(hax1,'Layer','top')
    set(gcf, 'PaperPositionMode', 'manual');
    print('-depsc',[dir_out prod '_Fig_1'],'-r300') % save figure 1
    close all
    
    % -------------------- plot latency histogram for month -------------------
    figure(2)
    xx = laten;
    if strcmp(data_type,'SIR_FDM_L2')
        if max(xx) < 9
            xx(abs(xx) > 8) = NaN;
            xhist = xhist1;
        else
            xx(abs(xx) > 14) = NaN;
            xhist = xhist2;
        end
    elseif strcmp(data_type,'SIR_IOP_L2')
        xx(abs(xx) > 4) = NaN;
    elseif strcmp(data_type,'SIR_GOP_L2')
        xx(abs(xx) > 36) = NaN;
    end
    %     got_to_this_point = 2 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point
    
    
    set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 7.4 6.3]);
    set(gca,'Units','centimeters','Position',[1.1 1.1 6 5])
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
else
    fprintf(fid,'%3.1f %3.1f %3.1f\n',[0 0 0]);
            fprintf(fid2,'%s\n','latency_fail');
        fprintf(fid2,'%3.1f\n',0);

end

% got_to_this_point = 3 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point

%% Data coverage and completeness
% total number of records for each day in last 30 days
nrec = data.nrec;

%clear t0f
fprintf(fid,'%-s\n',data.tStart{1});
fprintf(fid,'%-s\n',data.tStop{end});
fprintf(fid,'%u %u\n',[orbit0(1) orbitf(end)]);

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
perc_outSARinOrb = perc_outSARinOrb(io);
nth_orbit = 6324; %(5344 orbits/369 days and 1 record every 0.9434 s)
if strcmp(data_type,'SIR_FDM_L2')
    perc_orbit = nrecOrbit./(nth_orbit.*perc_LRMOrb./100).*100;
else
    perc_orbit = nrecOrbit./(nth_orbit.*perc_outSARinOrb./100).*100; % percentage records per orbit
end
perc_orbit = min(perc_orbit,100);
perc_orbitOcean = nrecOrbitOcean./(nth_orbit.*perc_orbOcean./100).*100;
perc_orbitOcean = min(perc_orbitOcean,100);

korb = find(orbitU == orbitUni(1));
orb_thr = 80;
% got_to_this_point = 2.99 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point


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
nth_ocean_noPolar = nth.*data.perc_ocean_noPolar./100;

num_rec = sum(nrec); % number of 1-Hz records for present day
if strcmp(data_type,'SIR_FDM_L2')
    nth = nth.*(perc_LRMDay./100);
    perc_rec = num_rec/sum(nth)*100;
else
    nth = nth.*(perc_outSARin./100);
    perc_rec = num_rec/sum(nth)*100; % precentage over theoretical max
end
perc_rec = min(perc_rec,100);
nthOcean = round(sum(nth_ocean));
nthTotal = round(sum(nth));
num_rec = min(num_rec,nthTotal);
fprintf(fid,'%u\n',nthTotal);
fprintf(fid,'%u\n',nthOcean);
fprintf(fid,'%u %4.1f\n',[num_rec perc_rec]);
% got_to_this_point = 2.66 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point


nrec_ocean = data.nrec_ocean;
perc_recOcean = nrec_ocean./nth_ocean.*100; %percentage over theoretical max
perc_recOcean = min(perc_recOcean,100);
nOceans = min(sum(nrec_ocean),nthOcean);
percOcean = min(nOceans/nthOcean*100,100);
fprintf(fid,'%u %4.1f\n',[nOceans percOcean]);

% generate warning if needed
ocean_thr = 80; %threshold
if perc_recOcean(end) < ocean_thr
    fprintf(fid2,'%s\n','ocean_dropout');
    fprintf(fid2,'%2.1f\n',perc_recOcean(end));
end

ylab = 'Percentage relative to theory (%) ';
leg = {'Total ','Ocean/lake '}; %#ok
pPos = [10.1 6.1];
pos = [8.7 5];

% plot percentage of records for each day in the month
pt = min(nrec./nth.*100,100);
mut = nanmean(pt);
stdt = nanstd(pt);
muo = nanmean(perc_recOcean);
stdo = nanstd(perc_recOcean);
leg1 = {['Total (\mu = ' sprintf('%0.1f',mut), ...
    ' , \sigma = ' sprintf('%0.1f',stdt) ')'], ...
    ['Ocean/lake (\mu = ' sprintf('%0.1f',muo), ...
    ' , \sigma = ' sprintf('%0.1f',stdo) ')']};
[hp1,hp2] = plot_timeSeries(2,(1:ndays)',pt,perc_recOcean,rgb1,rgb2,ylab, ...
    dates,leg1,pPos,pos);
delete([hp1 hp2])
rotateXLabels(gca(),45)
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out prod '_Fig_3'],'-r300') % save figure 2
close all
% got_to_this_point = 2.77 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point


% plot percentage of records per orbit in the month
mut = nanmean(perc_orbit);
stdt = nanstd(perc_orbit);
muo = nanmean(perc_orbitOcean);
stdo = nanstd(perc_orbitOcean);
leg1 = {['Total (\mu = ' sprintf('%0.1f',mut), ...
    ' , \sigma = ' sprintf('%0.1f',stdt) ')'], ...
    ['Ocean/lake (\mu = ' sprintf('%0.1f',muo), ...
    ' , \sigma = ' sprintf('%0.1f',stdo) ')']};
[hp1,hp2] = plot_timeSeries(3,orbitU,perc_orbit,perc_orbitOcean,rgb1, ...
    rgb2,ylab,orbitU,leg1,pPos,pos,'-');
delete([hp1 hp2])
%ylim([-15 100])
xlabel('Absolute orbit number','FontSize',8)
set(gca,'Layer','top')
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out prod '_Fig_4'],'-r300') % save figure 3
close all
% -------------------------------------- end data coverage and completeness
% got_to_this_point = 2.88 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point



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
% got_to_this_point = 4 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point
% got_to_this_point = 3 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point

cnt = 4; % counter for figures (three figure have already been saved)
for j=1:length(params) % loop over all parameters
    if j == 5
        j
    end
    %     got_to_this_point = 99 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point params j
    
    % ----------------- calculate percentage of edited by flag ----------------
    % this represents percentage of (total) records eddited over oceans and lakes,
    % and excluding polar regions.
    f1 = data.surface_type;
    f2 = data.(['validFlag_' vari{j}]);
    f1 = (f1 == 0 | f1 == 1);
    editFlag.(vari{j}) = (sum(f1 & data.noPolar) - sum(f2 & data.noPolar)) / ...
        sum(f1 & data.noPolar)*100;
    % -------------------------------------------------------------------------
    
    % ---------------------- Validity by NOC --------------------------
    nan_x = data.(['validFlag_' vari{j}]) & data.noPolar;
    x = params{j};
    ascFlag = logical(data.ascendingFlag); % used to plot ascending and descending
    if j ~= 5
        nan_x2 = data.(['validNOC_' vari{j}]);
        nrecx = data.(['nrecNOC_' vari{j}]); % total number of valid 1-Hz
    else
        nrecx = data.(['nrec_noPolar_' vari{j}]); % total number of valid 1-Hz
        nrecx2 = data.(['nrecPos_noPolar_' vari{j}]); % total number of valid 1-Hz
        nan_x2 = true(size(nan_x)) & x >= 0;
        percRecxTheorPos = (nrecx2./nth_ocean_noPolar).*100;
    end
    %         lonout = data.lon(~nan_x2 & nan_x);
    %         latout = data.lat(~nan_x2 & nan_x);
    %         recDateOut = recDate(~nan_x2 & nan_x);
    percRecxTheor = (nrecx./nth_ocean_noPolar).*100;
    nan_x2 = nan_x2 & nan_x;
    x = x(nan_x2);
    lonx = data.lon(nan_x2);
    latx = data.lat(nan_x2);
    depth = data.depth(nan_x2);
    ascFlag = ascFlag(nan_x2);
    %         recDatex = recDate(nan_x2);
    %         [recDatexU,ix] = unique(recDatex,'first');
    %         [~,ixout] = unique(recDateOut,'first');
    %         % save data
    %         fprintf(fid,'%u\n',num_recx);
    %         fprintf(fid,'%.1f\n',percRecx(end));
    %         fprintf(fid,'%.1f\n',percRecxTheor(end));
    
    
    % plot Percentage over theoretical max records for each day
    cnt = cnt+1;
    units = [{'SSH anomaly '};{'SWH '};{'Sigma0 '};{'Wind speed'};...
        {'off-nadir angle '}]; %#ok
    
    ylab = ['Percentage relative to theory (%) ']; %#ok
    pPos = [10 6.1];
    pos = [8.7 5];
    
    percRecxTheor(percRecxTheor > 100) = 100;
    
    mua = nanmean(percRecxTheor);
    stda = nanstd(percRecxTheor);
    if j == 5 && (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        percRecxTheorPos(percRecxTheorPos > 100) = 100;
        mub = nanmean(percRecxTheorPos);
        stdb = nanstd(percRecxTheorPos);
        leg1 = {['All values (\mu = ' sprintf('%0.1f',mua), ...
            ' , \sigma = ' sprintf('%0.1f',stda) ')'], ...
            ['Only positive values (\mu = ' sprintf('%0.1f',mub), ...
            ' , \sigma = ' sprintf('%0.1f',stdb) ')']};
        [hp1,hp2] = plot_timeSeries(cnt,(1:ndays)',percRecxTheor,percRecxTheorPos,rgb2, ...
            rgb1,ylab,dates,leg1,pPos,pos);
        delete([hp1 hp2])
    else
        leg1 = {['\mu = ' sprintf('%0.1f',mua), ...
            ' , \sigma = ' sprintf('%0.1f',stda)]};
        hp1 = plot_percData(cnt,(1:ndays)',percRecxTheor,rgb2,ylab,dates,leg1,pPos,pos);
        delete(hp1)
    end
    rotateXLabels(gca(),45)
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
    close all
    
    %     got_to_this_point = 99 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j
    
    
    %plot geographical distribution of valid data for present day
    if(j == 5)
        asc = ascFlag;
        ascTitle = ' (ascending passes)';
    else
        lat2 = -65:1:65;
        [muw,stdw] = weightMean(x,latx,lat2,1);
        asc = true(size(x));
        ascTitle = '';
    end
    cnt = cnt+1;
    units = [{'SSH anomaly (cm) '};{'SWH (m) '};{'Sigma0 (dB) '};...
        {'Altimeter wind speed (m/s)'};{'Square of off-nadir angle (10^{-2} deg^2) '}];
    titles = {'SSH anomaly','SWH','Sigma0','Altimeter wind speed','Square of off-nadir angle'};
    cytick = [5,1,1,2,1];
    scale = [100,1,1,1,100];
    
    strCol = {'p5','p25','median','p75','p95','mean','std'};
    
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        cmax = [20,8,14,20,5];
        cmin = [-20,0,7,0,0];
        strCol = [{'Mode'} strCol]; %#ok
        inLRM = data.inLRM;
        inLRM = inLRM(nan_x2);
        
        % calculate percentiles, mean and std
        prcOut = prctile(x(~inLRM & asc),[5 25 50 75 95]);
        mu1 = mean(x(~inLRM & asc));
        std1 = std(x(~inLRM & asc));
        
        dat_tab2 = {prcOut(1),prcOut(2),prcOut(3),prcOut(4) ...
            ,prcOut(5),mu1,std1};
        dat_tab2 = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab2,'unif',false);
    else
        mtmp = fix(median(x)*scale(j));
        cmax = [20,8,14,20,5];
        cmin = [-20,0,7,0,0];
        cmax(1) = cmax(1) + mtmp;
        cmin(1) = cmin(1) + mtmp;
        inLRM = true(size(x));
    end
    
    
    prcIn = prctile(x(inLRM & asc),[5 25 50 75 95]);
    mu1 = mean(x(inLRM & asc));
    std1 = std(x(inLRM & asc));
    
    % calculate geophysical global mean (no distinction LRM-PLRM)
    if j ~= 5
        strCol3 = {'Wmean','Wstd'};
        dat_tab3 = {muw,stdw};
        dat_tab3 = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab3,'unif',false);
    end
    dat_tab = {prcIn(1),prcIn(2),prcIn(3),prcIn(4) ...
        ,prcIn(5),mu1,std1};
    
    dat_tab = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab,'unif',false);
    titl = ['Science-valid ' titles{j} ' ' report_date ascTitle];
    pPos = [13.9 9.2];
    pos = [10.5 9];
    % % % % % % % % % % %     plot_map_nodisplay_monthly(cnt,lonx(asc),latx(asc),x(asc).*scale(j),cmax(j), ...
    % % % % % % % % % % %         cmin(j),cytick(j),units{j},titl,pPos,pos);
    % CHNAGED FOR TDS
    % figure(cnt)
    % fastscatterm(lonx(asc),latx(asc),x(asc).*scale(j),...
    %     [cmin(j) cmax(j)],[-180 180],[-90 90],3); title(titl)
    % load tempmonday.mat
    plot_map_nodisplay_monthly_v2(cnt,lonx(asc),latx(asc),x(asc).*scale(j),cmax(j), ...
        cmin(j),cytick(j),units{j},titl,pPos,pos);
    
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        xposTab = [4.5,4.5,4.5,4.5,4.5];
        sr = {'LRM','PLRM'};
        ht = plotTable(sr,strCol,[dat_tab;dat_tab2],[xposTab(j) 0.48],7,1);
        set(ht(:,1),'FontWeight','bold','BackGroundColor', ...
            [206/255 235/255 251/255])
        plotm(polar(2,:),polar(1,:),'color',[0 0 0],'LineWidth',1)
    else
        xposTab = [4.5,4.5,4.5,4.5,4.5];
        ht = plotTable([],strCol,dat_tab,[xposTab(j) 0.48],7,1);
    end
    %     got_to_this_point = 98 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j
    
    if j~=5
        ht2 = plotTable([],strCol3,dat_tab3,[2 0.48],7,1);
        set(ht2(1,:),'FontWeight','bold','BackGroundColor', ...
            [253/255 237/255 208/255])
    end
    set(ht(1,1:end-2),'FontWeight','bold','BackGroundColor', ...
        [206/255 235/255 251/255])
    set(ht(1,end-1:end),'FontWeight','bold','BackGroundColor', ...
        [253/255 237/255 208/255])
    
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
    close all
    
    % ----------- plot geographical distribution for descending node ----------
    if j == 5
        asc = ~ascFlag;
        ascTitle = ' (descending passes)';
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            
            % calculate percentiles, mean and std
            prcOut = prctile(x(~inLRM & asc),[5 25 50 75 95]);
            mu1 = mean(x(~inLRM & asc));
            std1 = std(x(~inLRM & asc));
            dat_tab2 = {prcOut(1),prcOut(2),prcOut(3),prcOut(4) ...
                ,prcOut(5),mu1,std1};
            dat_tab2 = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab2,'unif',false);
        end
        prcIn = prctile(x(inLRM & asc),[5 25 50 75 95]);
        mu1 = mean(x(inLRM & asc));
        std1 = std(x(inLRM & asc));
        dat_tab = {prcIn(1),prcIn(2),prcIn(3),prcIn(4) ...
            ,prcIn(5),mu1,std1};
        dat_tab = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab,'unif',false);
        titl = ['Science-valid ' titles{j} ' ' report_date ascTitle];
        pPos = [13.9 9.2];
        pos = [10.5 9];
        % CHANGED FOR TDS
        plot_map_nodisplay_monthly_v2(cnt,lonx(asc),latx(asc),x(asc).*scale(j),cmax(j), ...
            cmin(j),cytick(j),units{j},titl,pPos,pos);
        
        %         plot_map_nodisplay_monthly(cnt,lonx(asc),latx(asc),x(asc).*scale(j),cmax(j), ...
        %             cmin(j),cytick(j),units{j},titl,pPos,pos);
        %         if ~isempty(ixout)
        %             plotm(latout,lonout,'o','MarkerSize',2.5,'Color',[0.4 0.4 0.4])
        %         end
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            xposTab = [4.5,4.5,4.5,4.5,4.5];
            sr = {'LRM','PLRM'};
            ht = plotTable(sr,strCol,[dat_tab;dat_tab2],[xposTab(j) 0.48],7,1);
            set(ht(:,1),'FontWeight','bold','BackGroundColor', ...
                [206/255 235/255 251/255])
            plotm(polar(2,:),polar(1,:),'color',[0 0 0],'LineWidth',1)
        else
            ht = plotTable([],strCol,dat_tab,[xposTab(j) 0.48],8,1);
        end
        set(ht(1,1:end-2),'FontWeight','bold','BackGroundColor', ...
            [206/255 235/255 251/255])
        set(ht(1,end-1:end),'FontWeight','bold','BackGroundColor', ...
            [253/255 237/255 208/255])
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt) '_descend'],'-r300') % save figure
        close all
        
    end
    % -------------------------------------------------------------------------
    
    
    % ------------------------------- plot histogram --------------------------
    cnt = cnt+1;
    xtmp = params{j};
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        nan_x1 = nan_x2 & data.inLRM;
        nan_X2 = nan_x2 & ~data.inLRM;
        xx = xtmp(nan_x1).*scale(j);
        xx2 = xtmp(nan_X2).*scale(j);
        mu2 = mean(xx2);
        std2 = std(xx2);
    else
        xx = xtmp(nan_x2).*scale(j);
        xx2 = xx.*scale(j);
    end
    mu1 = mean(xx);
    std1 = std(xx);
    
    xlab = [{'SSH anomaly (cm)'};{'SWH (m)'};{'Sigma0 (dB)'}; ...
        {'Altimeter wind speed (m/s)'};{'Square of off-nadir angle (10^{-2} deg^2)'}];
    xmax = [50,12,20,24,6];
    
    xmin = [-50,0,4,0,0];
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        xmax(5) = 8;
    end
    xvalues = {-50:2:50,0:0.25:12,4:0.5:20,0:0.5:24,0:0.2:xmax(5)};
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
            xx(xx > 24 | xx == 0) = NaN;
            xx2(xx2 > 24 | xx2 == 0) = NaN;
        case 5
            xx(xx > xmax(5) | xx < xmin(j) | xx == 0) = NaN;
            xx2(xx2 > xmax(5) | xx2 < xmin(j) | xx2 == 0) = NaN;
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
    %     got_to_this_point = 97 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j
    
    px = get(gca,'XTick');
    py = get(gca,'YTick');
    Dx = px(end)-px(1);
    Dy = py(end)-py(1);
    units = [{'cm '};{'m '};{'dB '};{'m/s '};{'10^{-2} deg^2 '}];
    set(gca,'YTickLabel', py)
    
    if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        if j ~= 5
            fs = 3;
        else
            fs = 2.4;
        end
        text(px(end)-Dx/fs,py(end)-Dy/3.9,['\mu = ', ...
            sprintf('%0.1f',mu1) ' ' units{j}],'FontSize',8,'color',rgb2);
        text(px(end)-Dx/fs,py(end)-Dy*(1/3.9+0.07),['\sigma = ', ...
            sprintf('%0.1f',std1) ' ' units{j}],'FontSize',8,'color',rgb2);
        text(px(end)-Dx/fs,py(end)-Dy*(1/3.9+0.17),['\mu = ', ...
            sprintf('%0.1f',mu2) ' ' units{j}],'FontSize',8,'color',rgb1);
        text(px(end)-Dx/fs,py(end)-Dy*(1/3.9+0.24),['\sigma = ', ...
            sprintf('%0.1f',std2) ' ' units{j}],'FontSize',8,'color',rgb1);
    else
        if j ~= 5
            fs = 3.2;
        else
            fs = 2.4;
        end
        text(px(end)-Dx/fs,py(end)-Dy/12,['\mu = ', ...
            sprintf('%0.1f',mu1) ' ' units{j}],'FontSize',8);
        text(px(end)-Dx/fs,py(end)-Dy*(1/12+0.1),['\sigma = ', ...
            sprintf('%0.1f',std1) ' ' units{j}],'FontSize',8);
    end
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
    close all
    % -------------------------------------------------------------------------
    
    % ----------------plot geographical distribution of noise -----------------
    if j < 4
        scaleNoise = [100,100,100];
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            stdX1 = (data.(['noiseNOC_inLRM_' vari{j}])).*scaleNoise(j);
            %meanStdX1 = (data.(['noiseavgNOC_inLRM_' vari{j}])).*scaleNoise(j);
            stdX2 = (data.(['noiseNOC_outLRM_' vari{j}])).*scaleNoise(j);
            %meanStdX2 = (data.(['noiseavgNOC_outLRM_' vari{j}])).*scaleNoise(j);
            %         meanStdXNoc1 = (data.(['noiseavgNOC_inLRM_' vari{j}])).*scaleNoise(j);
            %         noiseStdNoc1 = (data.(['noiseStdNOC_inLRM_' vari{j}])).*scaleNoise(j);
            %         meanStdXNoc2 = (data.(['noiseavgNOC_outLRM_' vari{j}])).*scaleNoise(j);
            %         noiseStdNoc2 = (data.(['noiseStdNOC_outLRM_' vari{j}])).*scaleNoise(j);
        end
        stdX = (data.(['noiseNOC_' vari{j}])).*scaleNoise(j);
        %     meanStdX = (data.(['noiseavgNOC_' vari{j}])).*scaleNoise(j);
        %     %noiseStd = (data.(['noiseStd_' vari{j}])).*scaleNoise(j);
        %     noiseDiff = (data.(['noiseDiff_' vari{j}])).*scaleNoise(j);
        %noiseDiffNoc = (data.(['noiseDiffNOC_' vari{j}])).*scaleNoise(j);
        %     meanStdXNoc = (data.(['noiseavgNOC_' vari{j}])).*scaleNoise(j);
        %     noiseStdNoc = (data.(['noiseStdNOC_' vari{j}])).*scaleNoise(j);
        %
        %     % save data
        %     if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        %         fprintf(fid,'%.1f\n',meanStdX1(end));
        %         fprintf(fid,'%.1f\n',meanStdX2(end));
        %         fprintf(fid,'%.1f\n',meanStdX1(end)/sqrt(20));
        %         fprintf(fid,'%.1f\n',meanStdX2(end)/sqrt(20));
        %         fprintf(fid,'%.1f\n',noiseDiff(end));
        %         fprintf(fid,'%.1f\n',meanStdXNoc1(end));
        %         fprintf(fid,'%.1f\n',meanStdXNoc2(end));
        %         fprintf(fid,'%.1f\n',meanStdXNoc1(end)/sqrt(20));
        %         fprintf(fid,'%.1f\n',meanStdXNoc2(end)/sqrt(20));
        %         fprintf(fid,'%.1f\n',noiseDiffNoc(end));
        %     else
        %         fprintf(fid,'%.1f\n',meanStdX(end));
        %         fprintf(fid,'%.1f\n',meanStdX(end)/sqrt(20));
        %         fprintf(fid,'%.1f\n',noiseDiff(end));
        %         fprintf(fid,'%.1f\n',meanStdXNoc(end));
        %         fprintf(fid,'%.1f\n',meanStdXNoc(end)/sqrt(20));
        %         fprintf(fid,'%.1f\n',noiseDiffNoc(end));
        %     end
        %
        %     if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
        %         meanStdXNoc = meanStdXNoc1;
        %         noiseStdNoc = noiseStdNoc1;
        %     end
        
        % plot geographical distribution of 20-Hz SLA noise for present day
        cnt = cnt+1;
        units = [{'SSH anomaly noise (cm)'};{'SWH noise (cm)'}; ...
            {'Sigma0 noise (10^{-2} dB)'}];
        titles = {'Science-valid SSH anomaly noise','Science-valid SWH noise', ...
            'Science-valid sigma0 noise'};
        cytick = [2,10,2];
        zin = stdX;
        zin1 = zin;
        xin = lonx;
        yin = latx;
        
        strCol = {'p5','p25','median','p75','p95','mean','std'};
        
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            cmax = [17,90,20];
            cmin = [4,30,6];
            zin1 = stdX1;
            zin2 = stdX2;
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
        titl = [titles{j} ' ' report_date];
        pPos = [13.9 9.2];
        pos = [10.5 9];
        if size(zin,1) ~= size(xin,1) && size(zin,2) ~= size(xin,2)
            % CHANGED FOR TDS
            %             plot_map_nodisplay_monthly(cnt,xin,yin,zin,cmax(j),cmin(j),cytick(j), ...
            %                 units{j},titl,pPos,pos);
            plot_map_nodisplay_monthly_v2(cnt,xin,yin,zin,cmax(j),cmin(j),cytick(j), ...
                units{j},titl,pPos,pos);
            
            if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                sr = {'LRM','PLRM'};
                ht = plotTable(sr,strCol,[dat_tab;dat_tab2],[4.5 0.48],7,1);
                set(ht(:,1),'FontWeight','bold','BackGroundColor', ...
                    [206/255 235/255 251/255])
                plotm(polar(2,:),polar(1,:),'color',[0 0 0],'LineWidth',1)
            else
                ht = plotTable([],strCol,dat_tab,[4.5 0.48],8,1);
            end
            
            set(ht(1,1:end-2),'FontWeight','bold','BackGroundColor', ...
                [206/255 235/255 251/255])
            set(ht(1,end-1:end),'FontWeight','bold','BackGroundColor', ...
                [253/255 237/255 208/255])
            set(gcf, 'PaperPositionMode', 'manual')
            print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        end
        close all
    end
    % -------------------------------------------------------------------------
    
    % ---------------------- plot time series of noise ------------------------
    if j < 4
        cnt = cnt+1;
        units = [{'SSH anomaly noise (cm)'};{'SWH noise (cm)'}; ...
            {'Sigma0 noise (10^{-2} dB)'}];
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            meanStdXNoc = (data.(['noiseavgNOC_inLRM_' vari{j}])).*scaleNoise(j);
            noiseStdNoc = (data.(['noiseStdNOC_inLRM_' vari{j}])).*scaleNoise(j);
            meanStdXNoc2 = (data.(['noiseavgNOC_outLRM_' vari{j}])).*scaleNoise(j);
            noiseStdNoc2 = (data.(['noiseStdNOC_outLRM_' vari{j}])).*scaleNoise(j);
            ymax = max(max(meanStdXNoc),max(meanStdXNoc2)) ...
                +max(max(noiseStdNoc),max(noiseStdNoc2));
            ymin = min(min(meanStdXNoc),min(meanStdXNoc2)) ...
                -max(max(noiseStdNoc),max(noiseStdNoc2));
            dy = ymax-ymin;
            ylims = [ymin-0.1*dy ymax+0.1*dy];
        else
            meanStdXNoc = (data.(['noiseavgNOC_' vari{j}])).*scaleNoise(j);
            noiseStdNoc = (data.(['noiseStdNOC_' vari{j}])).*scaleNoise(j);
            ymax = max(meanStdXNoc)+max(noiseStdNoc);
            ymin = min(meanStdXNoc)-max(noiseStdNoc);
            dy = ymax-ymin;
            ylims = [ymin-0.1*dy ymax+0.1*dy];
        end
        ytick = {4:0.2:20,30:1:80,0:0.5:20};
        ytick2 = {4:2.0:20,25:10:160,0:2:26};
        hold on
        hs1 = errorbar((1:ndays)',meanStdXNoc,noiseStdNoc,'o','color',rgb2, ...
            'LineWidth',0.7,'MarkerSize',4,'MarkerFaceColor',rgb2);
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            if exist('xp','var')
                hold on
                hs2 = shadedErrorBar(xp,meanStdXNoc2,noiseStdNoc2, ...
                    {'-o','color',rgb1,'LineWidth', ...
                    1.5,'markerfacecolor',rgb1,'MarkerSize',3});
                set(hs2.edge,'LineStyle','-','color',rgb1,'LineWidth',0.5)
                set(hs2.patch,'Facecolor','none')
                hl = legend([hs1,hs2.mainLine],'LRM','PLRM','Location','northoutside',...
                    'Orientation','horizontal');
                set(hl, 'Box', 'off')
                set(hl, 'Color', 'none')
            end
        end
        ylabel(units{j},'FontSize',10)
        set(gca,'LineWidth',1,'FontSize',10,'XTick',unique([1:5:ndays ndays]), ...
            'XTickLabel',unique([dates(1:5:end);dates(end)]))
        xlim([0 ndays+1])
        box on
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 10.3 7.1]);
        %[0.5 0.5 10.3 6.7]
        set(gca,'Units','centimeters','Position',[1.4 1.25 8.7 5])
        if  ~isnan(mean(ylims))
            ylim(ylims)
        end
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
    %     got_to_this_point = 96 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j
    
    % -------------------------------------------------------------------------
    
    
    
    % --------plot performance as a function of distance from coast -----------
    if j < 4
        dmax = 30;
        Ylim1 = {[5,7],[40,50],[8.5,14]};
        Ylim2 = {[8.5,11],[50,80],[8.5,14]};
        dYTick = [0.5,5,1];
        binrange = (0.5:1:dmax+0.5)';
        xdist = (1:1:dmax)';
        d1 = data.distCoast20Hz./1000; % convert to km
        d1 = d1(:,nan_x2);
        e1 = data.(['absDiff_' vari{j}]).*100; % convert to cm or 10-2 dB
        kout = d1 > 1e15;
        d1(kout) = [];
        e1(kout) = [];
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            inLRM20 = data.inLRM20Hz;
            inLRM20 = inLRM20(:,nan_x2);
            inLRM20(kout) = [];
            d2 = d1(~inLRM20); %PLRM
            d1 = d1(inLRM20); %LRM
            e2 = e1(~inLRM20); %PLRM
            e1 = e1(inLRM20); %LRM
            
            [bcount,ind2] = histc(d2,binrange);
            E2 = NaN(1,length(bcount)-1); % skip the bin at the edge
            S2 = NaN(1,length(bcount)-1); % standard deviation for errobar
            for i=1:length(bcount)-1
                E2(i) = nanmedian(e2(ind2 == i));
                S2(i) = nanstd(e2(ind2 == i));
            end
        end
        [bcount,ind1] = histc(d1,binrange);
        E1 = NaN(1,length(bcount)-1);
        S1 = NaN(1,length(bcount)-1);
        for i=1:length(bcount)-1
            E1(i) = nanmedian(e1(ind1 == i));
            S1(i) = nanstd(e1(ind1 == i));
        end
        
        
        cnt = cnt+1;
        ylab = [{'20Hz SSH difference noise (cm) '};{'20Hz SWH difference noise (cm) '}; ...
            {'20Hz Sigma0 difference noise (10^{-2} dB) '}];
        pPos = [10 6.1];
        pos = [8.7 5];
        
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            plot_noiseDistCoast(cnt,xdist,E1,E2,rgb2,rgb1,ylab{j},pPos,pos,dmax, ...
                Ylim1{j},Ylim2{j},dYTick(j));
        else
            Ylim1 = {[5,6.5],[40,47],[8.5,12]};
            dYTick = [0.5,2,0.5];
            plot_noiseDistCoast(cnt,xdist,E1,[],rgb2,rgb1,ylab{j},pPos,pos,dmax, ...
                Ylim1{j},Ylim2{j},dYTick(j));
        end
        xlabel('Distance from coast (km)','FontSize',8)
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
        
    end
    % -------------------------------------------------------------------------
    
    
    
    % --------------------------- plot 2D histogram ---------------------------
    if j==1
        y_swh = data.swh;
        y_noise2DHist_inLRM = data.noise_ssha4scatter_inLRM;
        nan_2DHist_inLRM = data.nan2DHist_inLRM;
        f1 = data.validNOC_ssha & data.validNOC_swh;
        f3 = data.validFlag_ssha;
        f4 = data.validFlag_swh;
        f5 = data.noPolar;
        
        nan_2DHist_inLRM = nan_2DHist_inLRM & f1;
        if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            y_noise2DHist_outLRM = data.noise_ssha4scatter_outLRM;
            nan_2DHist_outLRM = data.nan2DHist_outLRM;
            nan_2DHist_outLRM = nan_2DHist_outLRM & f1;
            fin = f1(f3 & f4 & f5 & data.inLRM);
            fout = f1(f3 & f4 & f5 & ~data.inLRM);
            % CHANGED FOR TDS             y_noise2DHist_outLRM = y_noise2DHist_outLRM(fout);
        else
            fin = f1(f3 & f4 & f5);
        end
        % CHANGED FOR TDS    y_noise2DHist_inLRM = y_noise2DHist_inLRM(fin);
        
        
        % plot 2D histogram of noise in SLA vs SWH
        y = y_swh; % values for present day only
        nan_xy = nan_2DHist_inLRM; % values for present day only
        y = y(nan_xy);
        ssh_noise = y_noise2DHist_inLRM; % values for present day only
        ssh_noise(y == 0) = NaN;
        y(y == 0) = NaN;
        cnt = cnt+1;
        figure(cnt);
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 11.5 9.5]);
        set(gca,'Units','centimeters','Position',[1.2 1.1 10 8])
        axis_hist = [0 10 0 20];
        % CHANGED FOR TDS    z = [y' ssh_noise'.*100];
        % CHANGED FOR TDS     [cx,xmedian] = hist2d(z,axis_hist,133,267);
        hold on
        % CHANGED FOR TDS        plot(cx{1},xmedian,'LineWidth',2,'color','black')
        box on
        pu = get(gcf,'PaperUnits');
        pp = get(gcf,'PaperPosition');
        set(gcf,'Units',pu,'Position',pp)
        set(gca,'LineWidth',1.5,'FontSize',10)
        ylabel('SSH anomaly noise (cm)')
        xlabel('SWH (m)')
        xlim([0 10]) %30-day window
        %xlim([0 8]) % only 1 day
        ylim([0 20]) % 30-day window
        %ylim([1 18])
        set(gca,'XTick',0:2:16,'YTick',0:2:20)
        hc = colorbar;
        xpos = get(gca,'position');
        set(hc,'YTick',0:100:5000) % 30-day window
        %set(hc,'YTick',0:10:100) % only 1 day
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
            % CHANGED FOR TDS       z = [y' ssh_noise'.*100];
            z = [y' y'+NaN ] ;
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
            xlim([0 10]) %30-day window
            %xlim([0 8]) % only 1 day
            ylim([0 20]) % 30-day window
            %ylim([3 20])
            set(gca,'XTick',0:2:16,'YTick',0:2:20)
            hc = colorbar;
            xpos = get(gca,'position');
            set(hc,'YTick',0:20:1000) % 30-day window
            %set(hc,'YTick',0:2:20) % only 1 day
            set(hc,'FontSize',10)
            cpos = get(hc,'Position');
            cpos(3) = 0.6.*cpos(3);
            cpos(2) = cpos(2)+0.708;
            set(hc,'Position',cpos)
            set(gca,'position',xpos)
            text(6.7,4,'PLRM')
            set(gcf, 'PaperPositionMode', 'manual')
            print('-depsc',[dir_out prod '_Fig_' num2str(cnt) 'b'],'-r300') % save figure
            close all
        end
    end
    % -------------------------------------------------------------------------
    %     got_to_this_point = 95 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j
    
    
    % ---------------plot scatter at the edge of Pacific SAR box --------------
    if (~isempty(intersect(j,[1,2,4])) && (strcmp(data_type,'SIR_IOP_L2') ||...
            strcmp(data_type,'SIR_GOP_L2')))
        Scale = [100,1,0,1];
        Scale2 = [1,100,0,100];
        offSet = [2,0.2,0,1];
        kin = data.kInPacificSAR;
        kout = data.kOutPacificSAR;
        f = data.(['validFlag_' vari{j}]) & data.(['validNOC_' vari{j}]);
        f = f(kin) & f(kout);
        kin = kin(f);
        kout = kout(f);
        dks = kout - kin;
        krms = kout + dks;
        x = params{j};
        TTi = data.tutc.*24*3600;
        Ti = TTi(kin);
        To = TTi(kout);
        Tr = TTi(krms);
        xIn = x(kin).*Scale(j);
        xOut = x(kout).*Scale(j);
        xRms = x(krms).*Scale(j);
        Xmin = min(min(xIn),min(xOut));
        Xmax = max(max(xIn),max(xOut));
        Xmin = max(Xmin,-40);
        Xmax = min(Xmax,40);
        
        si = nanstd(xIn - xOut);
        sr = nanstd(xOut - xRms);
        ki = (abs(xIn - xOut) < 3.*si) & (abs(Ti - To) < 1);
        kr = (abs(xOut - xRms) < 3.*sr) & (abs(To - Tr) < 1);
        di = (xIn - xOut).^2;
        dr = (xOut - xRms).^2;
        
        % compute RMS of differences
        rmsEdge = sqrt(nanmean(di(ki)));
        rmsOut = sqrt(nanmean(dr(kr)));
        
        fprintf(fid,'%4.1f\n',rmsEdge.*Scale2(j));
        fprintf(fid,'%4.1f\n',rmsOut.*Scale2(j));
        
        cnt = cnt+1;
        unitsX = [{'SSH anomaly LRM (cm)'};{'SWH LRM (m)'};{' '}; ...
            {'Wind speed LRM (m/s)'}];
        unitsY = [{'SSH anomaly PLRM (cm)'};{'SWH PLRM (m)'};{' '}; ...
            {'Wind speed PLRM (m/s)'}];
        Xlim = [Xmin-offSet(j),Xmax+offSet(j)];
        xTick = {-50:10:50,0:1:30,[],0:5:30};
        
        pPos = [10 10];
        pos = [8 8];
        figure(cnt);
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
        set(gca,'Units','centimeters','Position',[1.2 1 pos])
        scatter(xOut,xIn,40,rgb1,'filled')
        ylabel(gca,unitsY{j},'FontSize',8)
        xlabel(gca,unitsX{j},'FontSize',8)
        set(gca,'LineWidth',1,'FontSize',8,'XTick',xTick{j},'YTick',xTick{j})
        if ~isempty(Xlim)
            xlim(Xlim)
            ylim(Xlim)
        end
        hold on
        plot(-1e2:1e2,-1e2:1e2,'color','black','LineWidth',1.5)
        box on
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
    end
    % -------------------------------------------------------------------------
    
    
    
    %% ------------------ calculate crossover for data in month ---------------
    
    
    date0 = datenum(report_date,'mmm-yyyy');
    [Y,M] = datevec(date0);
    datef = date0+eomday(Y,M)-1;
    
    if j == 1 && ~exist([path3 'xOver_' data_type '_' report_date '.mat'],'file')
        maxd = 10; % time window (days) for the crossovers
        
        % Ensure data are available in matlab structure, if not process them
        dleft = date0-maxd:date0-1;
        dright = datef+1:datef+maxd;
        
        % -----------------------first process left data if available
        t_window = datestr(dleft,'yyyymmdd');
        ndays1 = size(t_window,1); % number of days in month
        files = [repmat([data_type '_'],ndays1,1) t_window repmat('.mat',ndays1,1)];
        kempty = false(1,size(files,1));
        for i=1:size(files,1)
            if ~exist([path2 files(i,:)],'file')
                read_cryo2Data(in_dir,path2,t_window(i,:),t_window(i,:));
            end
            if ~exist([path2 files(i,:)],'file')
                kempty(i) = true;
            end
        end
        
        % Read statistics for month
        files = [repmat(['Stats_' data_type '_'],ndays1,1) t_window repmat('.mat',ndays1,1)];
        dataLeft = struct([]);
        if any(~kempty)
            files = files(~kempty,:);
            for i=1:size(files,1)
                if ~exist([path1 files(i,:)],'file')
                    daily_stats_v6_TDS(t_window(i,:),data_type,path2,in_dir,path1, ...
                        DTU10MSS,gshhs_i,A);
                end
                dataLeft = [dataLeft load([path1 files(i,:)])]; %#ok
            end
        end
        ndays1 = sum(~kempty);
        % ------------------------ end processing left data
        
        
        % -----------------------now process right data if available
        t_window = datestr(dright,'yyyymmdd');
        ndays2 = size(t_window,1); % number of days in month
        files = [repmat([data_type '_'],ndays2,1) t_window repmat('.mat',ndays2,1)];
        kempty = false(1,size(files,1));
        for i=1:size(files,1)
            if ~exist([path2 files(i,:)],'file')
                read_cryo2Data(in_dir,path2,t_window(i,:),t_window(i,:));
            end
            if ~exist([path2 files(i,:)],'file')
                kempty(i) = true;
            end
        end
        
        % Read statistics for month
        files = [repmat(['Stats_' data_type '_'],ndays2,1) t_window repmat('.mat',ndays2,1)];
        dataRight = struct([]);
        if any(~kempty)
            files = files(~kempty,:);
            for i=1:size(files,1)
                if ~exist([path1 files(i,:)],'file')
                    daily_stats_v6_TDS(t_window(i,:),data_type,path2,in_dir,path1, ...
                        DTU10MSS,gshhs_i,A);
                end
                dataRight = [dataRight load([path1 files(i,:)])]; %#ok
            end
        end
        ndays2 = sum(~kempty);
        % ------------------------ end processing right data
        
        if ~isempty(dataLeft)
            dataLeft = rmfield(dataLeft,'perc_orbOcean');
            dataLeft = rmfield(dataLeft,'perc_LRMOrb');
            dataLeft = rmfield(dataLeft,'polarMask');
            dataLeft = rmfield(dataLeft,'perc_outSARinOrb');
            dataLeft = rmfield(dataLeft,'modeSeas');
            flds = fieldnames(dataLeft);
            dataLeft = cellfun(@(x) {horzcat(dataLeft.(x))},flds);
            dataLeft = cell2struct(dataLeft,flds);
            
            nan_xLeft = dataLeft.(['validFlag_' vari{j}]) & dataLeft.noPolar;
            nan_x2Left = dataLeft.(['validNOC_' vari{j}]);
            nan_x2Left = nan_xLeft & nan_x2Left;
            if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                inLRMLeft = dataLeft.inLRM;
                inLRMLeft = inLRMLeft(nan_x2Left);
            end
            tutcLeft = dataLeft.tutc(nan_x2Left);
            orbLeft = dataLeft.orbits(nan_x2Left);
            lonLeft = dataLeft.lon(nan_x2Left);
            latLeft = dataLeft.lat(nan_x2Left);
            nRLeft = dataLeft.(['nrecNOC_' vari{j}]);
            xLeft = dataLeft.(vari{j});
            xLeft = xLeft(nan_x2Left);
            depthLeft = dataLeft.depth(nan_x2Left);
        end
        
        if ~isempty(dataRight)
            dataRight = rmfield(dataRight,'perc_orbOcean');
            dataRight = rmfield(dataRight,'perc_LRMOrb');
            dataRight = rmfield(dataRight,'polarMask');
            dataRight = rmfield(dataRight,'perc_outSARinOrb');
            dataRight = rmfield(dataRight,'modeSeas');
            flds = fieldnames(dataRight);
            dataRight = cellfun(@(x) {horzcat(dataRight.(x))},flds);
            dataRight = cell2struct(dataRight,flds);
            
            nan_xRight = dataRight.(['validFlag_' vari{j}]) & dataRight.noPolar;
            nan_x2Right = dataRight.(['validNOC_' vari{j}]);
            nan_x2Right = nan_xRight & nan_x2Right;
            if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                inLRMRight = dataRight.inLRM;
                inLRMRight = inLRMRight(nan_x2Right);
            end
            tutcRight = dataRight.tutc(nan_x2Right);
            orbRight = dataRight.orbits(nan_x2Right);
            lonRight = dataRight.lon(nan_x2Right);
            latRight = dataRight.lat(nan_x2Right);
            nRRight = dataRight.(['nrecNOC_' vari{j}]);
            xRight = dataRight.(vari{j});
            xRight = xRight(nan_x2Right);
            depthRight = dataRight.depth(nan_x2Right);
        end
        
        if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            inLRMCross = data.inLRM(data.validFlag_ssha & data.validNOC_ssha & ...
                data.noPolar);
            if ndays1 ~= 0
                inLRMCross = [inLRMLeft inLRMCross]; %#ok
            end
            
            if ndays2 ~= 0
                inLRMCross = [inLRMCross inLRMRight]; %#ok
            end
            
        end
        
        tutc = data.tutc(nan_x2);
        orb = data.orbits(nan_x2);
        nR = data.(['nrecNOC_' vari{j}]);
        Lonx = lonx;
        Latx = latx;
        x = params{j};
        Xcross = x(nan_x2);
        if ndays1 ~= 0
            tutc = [tutcLeft tutc]; %#ok
            orb = [orbLeft orb]; %#ok
            nR = [nRLeft nR];%#ok
            Lonx = [lonLeft Lonx];%#ok
            Latx = [latLeft Latx];%#ok
            Xcross = [xLeft Xcross]; %#ok
            depth = [depthLeft depth]; %#ok
        end
        if ndays2 ~= 0
            tutc = [tutc tutcRight];%#ok
            orb = [orb orbRight];%#ok
            nR = [nR nRRight];%#ok
            Lonx = [Lonx lonRight];%#ok
            Latx = [Latx latRight];%#ok
            Xcross = [Xcross xRight]; %#ok
            depth = [depth depthRight]; %#ok
        end
        ndaysT = ndays+ndays1+ndays2;
        
        nRc = cumsum(nR);
        nRc = [1 nRc]; %#ok
        kcross = zeros(2,50000);
        tstep = maxd+1;
        iloop = unique([(maxd+1):tstep:(ndaysT-maxd) ndaysT-maxd]);
        % the iloop will avoid testing points for crossovers that have already
        % been tested
        for ii=iloop
            disp(ii)
            d0 = ii-maxd;
            df = ii+maxd;
            kday = nRc(d0):nRc(df+1);
            if kday ~= 0
                xin = Lonx(kday);
                yin = Latx(kday);
                oin = orb(kday);
                kx = crossOver(xin,yin,oin);
            end
            if exist('kx','var')
                if ~isempty(kx)
                    kx(1,:) = kday(kx(1,:));
                    kx(2,:) = kday(kx(2,:));
                    klast = find(kcross(1,:)==0,1,'first');
                    kcross(:,klast:klast+size(kx,2)-1) = kx;
                end
            end
        end
        klast = find(kcross(1,:)==0,1,'first');
        kcross = kcross(:,1:klast-1);
        
        % reject duplicated crossovers in a fast way
        [ksort,Isort] = sort(kcross,1);
        [kcross,Ia,~] = unique(ksort','rows');
        kchanged = find(Isort(1,:) == 2);
        [~,i1,~] = intersect(Ia,kchanged);
        kcross(i1,[1 2]) = kcross(i1,[2 1]); %reposition pairs
        kcross = kcross';
        
        % reject crossovers for points outside the month
        tout = (tutc(kcross(1,:)) < date0 & tutc(kcross(2,:)) < date0) | ...
            (tutc(kcross(1,:)) > datef+1 & tutc(kcross(2,:)) > datef+1);
        kcross(:,tout) = [];
        
        
        if all(size(kcross))
            lonCross = Lonx(kcross(1,:));
            latCross = Latx(kcross(1,:));
            tutcCross = tutc(kcross(1,:));
            if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                inLRMCross = inLRMCross(kcross(1,:));
            end
            depth = depth(kcross(1,:));
            diffCross = zeros(1,size(kcross,2));
            for z=1:size(kcross,2)
                dt1 = abs(tutc-tutc(kcross(1,z))).*(24*3600) <= 2;
                dt2 = abs(tutc-tutc(kcross(2,z))).*(24*3600) <= 2;
                dtx = abs(tutc(kcross(1,z)) - tutc(kcross(2,z)));
                if dtx < maxd % days
                    diffCross(1,z) = median(Xcross(dt1))-median(Xcross(dt2));
                else
                    continue
                end
            end
            lonCross(diffCross == 0) = [];
            latCross(diffCross == 0) = [];
            tutcCross(diffCross == 0) = [];
            if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                inLRMCross(diffCross == 0) = [];
            end
            depth(diffCross == 0) =[];
            diffCross(diffCross == 0) = [];
            
        else
            diffCross = NaN;
            lonCross = NaN;
            latCross = NaN;
            tutcCross = NaN;
            depth = NaN;
            if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                inLRMCross = NaN;
            end
        end
        xOver.diffCross = diffCross;
        xOver.lonCross = lonCross;
        xOver.latCross = latCross;
        xOver.tutcCross = tutcCross;
        if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            xOver.inLRMCross = inLRMCross;
        end
        xOver.depth = depth;
        save([path3 'xOver_' data_type '_' report_date '.mat'],'-struct','xOver')
    else
        load([path3 'xOver_' data_type '_' report_date '.mat'])
    end
    % -------------------------------------------------------------------------
    
    
    % ------------------------ plot crossover map -----------------------------
    if j == 1
        xin = lonCross;
        yin = latCross;
        zin = diffCross.*100;
        
        cnt = cnt+1;
        units = {'Crossover differences (cm)'};
        titles = {'Science-valid crossover differences'};
        
        strCol = {'p5','p25','median','p75','p95','mean','std'};
        cmax = 15;
        cmin = -15;
        cytick = 5;
        
        dat_tab = {prctile(zin,5),prctile(zin,25),prctile(zin,50), ...
            prctile(zin,75),prctile(zin,95),nanmean(zin),nanstd(zin)};
        dat_tab = cellfun(@(x) sprintf('%0.1f',x),dat_tab,'unif',false);
        titl = [titles{j} ' ' report_date];
        pPos = [13.9 9.2];
        pos = [10.5 9];
        % CHANGED FOR TDS
        plot_map_nodisplay_monthly_v2(cnt,xin,yin,zin,cmax,cmin,cytick, ...
            units{j},titl,pPos,pos,true);
        
        %         plot_map_nodisplay_monthly(cnt,xin,yin,zin,cmax,cmin,cytick, ...
        %             units{j},titl,pPos,pos,true);
        ht = plotTable([],strCol,dat_tab,[4 0.48],8,1);
        set(ht(1,1:end-2),'FontWeight','bold','BackGroundColor', ...
            [206/255 235/255 251/255])
        set(ht(1,end-1:end),'FontWeight','bold','BackGroundColor', ...
            [253/255 237/255 208/255])
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
    end
    %     got_to_this_point = 90 ;
    %     save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j
    
    % -------------------------------------------------------------------------
    
    
    % ------------------ plot time series of std of crossovers ----------------
    if j == 1
        tin = floor(tutcCross);
        % perform screening
        k20 = abs(zin) < 20 & abs(yin) <=50; % reject crossover differences larger than 20 cm and points > 50�lat
        zin_filt = zin(k20);
        tin_filt = tin(k20);
        depth_filt = depth(k20);
        if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            inLRM_filt = inLRMCross(k20);
        end
        
        k1000 = depth_filt < -1000; % reject points over shallow water
        zin_filt = zin_filt(k1000);
        tin_filt = tin_filt(k1000);
        if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            inLRM_filt = inLRM_filt(k1000);
        end
        
        % compute statistics
        days = date0:datef;
        std1 = NaN(1,ndays);
        std2 = NaN(1,ndays);
        if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            stdIn = NaN(1,ndays);
            stdOut = NaN(1,ndays);
        end
        for nn=1:ndays
            kday = (tin == days(nn));
            std1(nn) = nanstd(zin(kday));
            kday = (tin_filt == days(nn));
            std2(nn) = nanstd(zin_filt(kday));
            if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
                stdIn(nn) = nanstd(zin_filt(kday & inLRM_filt));
                stdOut(nn) = nanstd(zin_filt(kday & ~inLRM_filt));
            end
        end
        
        
        % ---------------------plot std of crossovers
        ylab = 'Standard deviation (cm) ';
        pPos = [10.1 6.1];
        pos = [8.7 5];
        
        % plot percentage of records for each day in the month
        cnt = cnt+1;
        if (strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_L2'))
            mut = nanmean(stdIn);
            stdt = nanstd(stdIn);
            muo = nanmean(stdOut);
            stdo = nanstd(stdOut);
            leg1 = {['LRM (\mu = ' sprintf('%0.1f',mut), ...
                ' , \sigma = ' sprintf('%0.1f',stdt) ')'], ...
                ['PLRM (\mu = ' sprintf('%0.1f',muo), ...
                ' , \sigma = ' sprintf('%0.1f',stdo) ')']};
            if ~isnan(nanmean(stdIn))
                [hp1,hp2] = plot_stdCrossover(cnt,(1:ndays)',stdIn,stdOut,rgb1,rgb2,ylab, ...
                    dates,leg1,pPos,pos);
            end
        else
            mut = nanmean(std2);
            stdt = nanstd(std2);
            leg1 = {['\mu = ' sprintf('%0.1f',mut), ...
                ' , \sigma = ' sprintf('%0.1f',stdt)]};
            if ~isnan(nanmean(stdIn))
                [hp1,hp2] = plot_stdCrossover(cnt,(1:ndays)',std2,[],rgb1,rgb2,ylab, ...
                    dates,[],pPos,pos);
                yl = ylim;
                text(1.3924,yl(2)-0.25,leg1,'FontSize',8);
            end
        end
        
        if ~isnan(nanmean(stdIn))
            
            delete([hp1 hp2])
        end
        rotateXLabels(gca(),45)
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300')
        close all
    end
    
    % -------------------------------------------------------------------------
    
end
% % 
% % % REPLACE fig 6 with version run from Matlab direct as for some weird
% % % reason this does not work in nodisplay mode
% % pwd_use = pwd ;
% % cd /noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/TDS
% % delete GOP_Fig_6.*
% % copyfile fig_6/GOP_Fig_6.eps 
% % cd(pwd_use) ; clear pwd_use


%% --------------------- plot corrections for GOP --------------------------
if strcmp(data_type,'SIR_GOP_L2')
    params = {data.wet,data.dry,data.iono,data.ssb,data.atm};
    f1 = data.surface_type;
    f1 = (f1 == 0 | f1 == 1);
    
    for j=1:length(params) % loop over all corrections
        % --------------------- ascending first ---------------------------
        ascFlag = logical(data.ascendingFlag);
        ascTitle = ' (ascending passes)';
        x = params{j};
        nan_x = data.noPolar & x == x & f1 & ascFlag;
        x = x(nan_x);
        lonx = data.lon(nan_x);
        latx = data.lat(nan_x);
        strCol = {'p5','p25','median','p75','p95','mean','std'};
        prcIn = prctile(x,[5 25 50 75 95]);
        
        %plot geographical distribution
        cnt = cnt+1;
        units = [{'Wet tropospheric correction (cm) '}; ...
            {'Dry tropospheric correction (m) '};{'Ionospheric correction (cm) '};...
            {'Sea state bias (cm)'};{'Atmospheric correction (cm) '}];
        titles = {'Wet tropospheric correction','Dry tropospheric correction', ...
            'Ionospheric correction','Sea state bias','Atmospheric correction'};
        cytick = [5,0.02,3,5,5];
        scale = [100,1,100,100,100];
        xposTab = [4,4,4,4,4];
        cmax = [0 -2.2 -3 0 26];
        cmin = [-40 -2.34 -12 -20 -26];
        cmax(3) = fix((median(x) + std(x)*3)*scale(3));
        cmin(3) = fix((median(x) - std(x)*3)*scale(3));
        dtmp = cmax(3) - cmin(3);
        dtmp = dtmp/6;
        if(round(dtmp) == 0)
            cytick(3) = round(dtmp*10)/10;
        else
            cytick(3) = round(dtmp);
        end
        
        mu1 = mean(x);
        std1 = std(x);
        dat_tab = {prcIn(1),prcIn(2),prcIn(3),prcIn(4) ...
            ,prcIn(5),mu1,std1};
        if j ~= 2
            dat_tab = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab,'unif',false);
        else
            dat_tab = cellfun(@(x) sprintf('%0.2f',x*scale(j)),dat_tab,'unif',false);
        end
        titl = [titles{j} ' ' report_date ascTitle];
        pPos = [13.9 9.4];
        pos = [10.5 9];
        % CHANGED FOR TDS plot_map_nodisplay_monthly(cnt,lonx,latx,x.*scale(j),cmax(j), ...
        %             cmin(j),cytick(j),units{j},titl,pPos,pos);
        plot_map_nodisplay_monthly_v2(cnt,lonx,latx,x.*scale(j),cmax(j), ...
            cmin(j),cytick(j),units{j},titl,pPos,pos);
        
        ht = plotTable([],strCol,dat_tab,[xposTab(j) 0.48],8,1);
        set(ht(1,1:end-2),'FontWeight','bold','BackGroundColor', ...
            [206/255 235/255 251/255])
        set(ht(1,end-1:end),'FontWeight','bold','BackGroundColor', ...
            [253/255 237/255 208/255])
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt)],'-r300') % save figure
        close all
        
        
        % ------------------------ Now descending -------------------------
        ascTitle = ' (descending passes)';
        x = params{j};
        nan_x = data.noPolar & x == x & f1 & ~ascFlag;
        x = x(nan_x);
        lonx = data.lon(nan_x);
        latx = data.lat(nan_x);
        
        %plot geographical distribution
        prcIn = prctile(x,[5 25 50 75 95]);
        cmax(3) = fix((median(x) + std(x)*3)*scale(3));
        cmin(3) = fix((median(x) - std(x)*3)*scale(3));
        dtmp = cmax(3) - cmin(3);
        dtmp = dtmp/6;
        if(round(dtmp) == 0)
            cytick(3) = round(dtmp*10)/10;
        else
            cytick(3) = round(dtmp);
        end
        
        mu1 = mean(x);
        std1 = std(x);
        dat_tab = {prcIn(1),prcIn(2),prcIn(3),prcIn(4) ...
            ,prcIn(5),mu1,std1};
        if j ~= 2
            dat_tab = cellfun(@(x) sprintf('%0.1f',x*scale(j)),dat_tab,'unif',false);
        else
            dat_tab = cellfun(@(x) sprintf('%0.2f',x*scale(j)),dat_tab,'unif',false);
        end
        titl = [titles{j} ' ' report_date ascTitle];
        pPos = [13.9 9.4];
        pos = [10.5 9];
        % CHANGED FOR TDS plot_map_nodisplay_monthly(cnt,lonx,latx,x.*scale(j),cmax(j), ...
        %             cmin(j),cytick(j),units{j},titl,pPos,pos);
        plot_map_nodisplay_monthly_v2(cnt,lonx,latx,x.*scale(j),cmax(j), ...
            cmin(j),cytick(j),units{j},titl,pPos,pos);
        
        ht = plotTable([],strCol,dat_tab,[xposTab(j) 0.48],8,1);
        set(ht(1,1:end-2),'FontWeight','bold','BackGroundColor', ...
            [206/255 235/255 251/255])
        set(ht(1,end-1:end),'FontWeight','bold','BackGroundColor', ...
            [253/255 237/255 208/255])
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' num2str(cnt) '_descend'],'-r300') % save figure
        close all
        
        
    end
    % -------------------------------------------------------------------------
    
    % ---------------------------- compute spectrum ---------------------------
    w = 128;
    w_hi = w*20;
    nd = length(spectra.klrm);
    P_lrm = cell(nd,1);
    F_lrm = cell(nd,1);
    P_plrm = cell(nd,1);
    F_plrm = cell(nd,1);
    for j=1:nd
        
        klrm = spectra.klrm{j};
        h_hi = spectra.h_hi{j};
        t_hi = spectra.t_hi{j};
        f1 = spectra.f1{j};
        f2 = spectra.f2{j};
        
        h_lrm = h_hi(:,f1 & f2 & klrm);
        t_lrm = t_hi(:,f1 & f2 & klrm);
        h_lrm = h_lrm(:);
        t_lrm = t_lrm(:);
        k = h_lrm == h_lrm;
        h_lrm = h_lrm(k);
        t_lrm = t_lrm(k);
        d_lrm = diff(t_lrm);
        
        h_plrm = h_hi(:,f1 & f2 & ~klrm);
        t_plrm = t_hi(:,f1 & f2 & ~klrm);
        h_plrm = h_plrm(:);
        t_plrm = t_plrm(:);
        k = h_plrm == h_plrm;
        h_plrm = h_plrm(k);
        t_plrm = t_plrm(k);
        d_plrm = diff(t_plrm);
        
        
        % for the 20Hz data most gaps are of 20 records
        k = find(d_lrm > 0.05); % median time separation is 0.0472 s
        k2 = find(d_lrm(k) < 2.1);
        ct = 0;
        for i=1:length(k2)
            ki = k(k2(i)) + ct;
            np = round(d_lrm(ki-ct)/0.0472) - 1;
            if(sum(diff(t_lrm) < 0) == 0)
                hi = interp1(t_lrm,h_lrm,t_lrm(ki)+(0.0472.*(1:np)'),'linear') +...
                    randn(np,1)*0.062; % add white noise
                t_lrm = [t_lrm(1:ki) ; t_lrm(ki)+(0.0472.*(1:np)') ; t_lrm(ki+1:end)];
                h_lrm = [h_lrm(1:ki) ; hi ; h_lrm(ki+1:end)];
            else
                np = 0;
            end
            ct = ct + np;
        end
        
        k = find(d_plrm > 0.05); % median time separation is 0.0472 s
        k2 = find(d_plrm(k) < 2.1);
        ct = 0;
        for i=1:length(k2)
            ki = k(k2(i)) + ct;
            np = round(d_plrm(ki-ct)/0.0472) - 1;
            if(sum(diff(t_plrm) < 0) == 0)
                hi = interp1(t_plrm,h_plrm,t_plrm(ki)+(0.0472.*(1:np)'),'linear') +...
                    randn(np,1)*0.062; % add white noise
                t_plrm = [t_plrm(1:ki) ; t_plrm(ki)+(0.0472.*(1:np)') ; t_plrm(ki+1:end)];
                h_plrm = [h_plrm(1:ki) ; hi ; h_plrm(ki+1:end)];
            else
                np = 0;
            end
            ct = ct + np;
        end
        
        
        n_lrm = length(h_lrm);
        d_lrm = diff(t_lrm);
        k_lrm = find(d_lrm > 0.05); % median time separarion is 0.0472
        k_lrm = [k_lrm; n_lrm]; %#ok
        k_lrm(2:end) = k_lrm(2:end) - k_lrm(1:end-1);
        h_lrm = mat2cell(h_lrm,k_lrm,1); % each cell contains data for a valid segment
        
        
        n_plrm = length(h_plrm);
        d_plrm = diff(t_plrm);
        k_plrm = find(d_plrm > 0.05); % median time separarion is 0.0472
        k_plrm = [k_plrm; n_plrm]; %#ok
        k_plrm(2:end) = k_plrm(2:end) - k_plrm(1:end-1);
        h_plrm = mat2cell(h_plrm,k_plrm,1); % each cell contains data for a valid segment
        
        np_lrm = cellfun(@length,h_lrm);
        h_lrm = h_lrm(np_lrm >= w_hi);
        n_lrm = length(h_lrm);
        
        np_plrm = cellfun(@length,h_plrm);
        h_plrm = h_plrm(np_plrm >= w_hi);
        n_plrm = length(h_plrm);
        
        pxx_lrm = cell(size(h_lrm));
        f_lrm = cell(size(h_lrm));
        for i=1:n_lrm
            [pxx_lrm{i},f_lrm{i}] = pwelch(h_lrm{i},w_hi,w_hi/2,w_hi,1/0.32);
        end
        np_lrm = cellfun(@length,f_lrm);
        f_lrm = cell2mat(f_lrm);
        pxx_lrm = cell2mat(pxx_lrm);
        
        pxx_plrm = cell(size(h_plrm));
        f_plrm = cell(size(h_plrm));
        for i=1:n_plrm
            [pxx_plrm{i},f_plrm{i}] = pwelch(h_plrm{i},w_hi,w_hi/2,w_hi,1/0.32);
        end
        np_plrm = cellfun(@length,f_plrm);
        f_plrm = cell2mat(f_plrm);
        pxx_plrm = cell2mat(pxx_plrm);
        
        if ~isempty(np_lrm)
            f_lrm = reshape(f_lrm,np_lrm(1),length(f_lrm)/np_lrm(1));
            pxx_lrm = reshape(pxx_lrm,np_lrm(1),length(pxx_lrm)/np_lrm(1));
        else
            f_lrm = [];
            pxx_lrm = [];
        end
        
        
        if ~isempty(np_plrm)
            f_plrm = reshape(f_plrm,np_plrm(1),length(f_plrm)/np_plrm(1));
            pxx_plrm = reshape(pxx_plrm,np_plrm(1),length(pxx_plrm)/np_plrm(1));
        else
            f_plrm = [];
            pxx_plrm = [];
        end
        
        P_lrm{j} = pxx_lrm;
        F_lrm{j} = f_lrm;
        P_plrm{j} = pxx_plrm;
        F_plrm{j} = f_plrm;
    end
    
    P_lrm = cell2mat(P_lrm');
    F_lrm = cell2mat(F_lrm');
    P_lrm = mean(P_lrm,2);
    F_lrm = mean(F_lrm,2);
    
    P_plrm = cell2mat(P_plrm');
    F_plrm = cell2mat(F_plrm');
    P_plrm = mean(P_plrm,2);
    F_plrm = mean(F_plrm,2);
    
    
    h1 = loglog(F_plrm(1:end-1),P_plrm(1:end-1),'color',rgb2,'linewidth',2);
    hold on
    h2 = loglog(F_lrm(1:end-1),P_lrm(1:end-1),'color',rgb1,'linewidth',2);
    save temp_debug.mat h1 h2
    if ~isempty(h1) && ~isempty(h2)
        legend([h2,h1],'20Hz CryoSat-2 LRM','20Hz CryoSat-2 PLRM')
    end
    set(gca,'FontSize',16)
    ylabel('PSD (m^2/cpkm)')
    xlabel('Wavenumber (cpkm)')
    grid on
    
    set(gcf, 'PaperPositionMode', 'manual')
    print('-depsc',[dir_out prod '_Fig_' num2str(cnt+1)],'-r300') % save figure
    close all
end
% -------------------------------------------------------------------------
% got_to_this_point = 1000000 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point

% ------------------ compute percentage of editted data -------------------
n1 = data.nrec_noPolar_ssha;
n2 = data.nrec_noPolar_swh;
n3 = data.nrec_noPolar_sigma0;
n4 = data.nrec_noPolar_wsp;
editSLA = sum(n1.*data.percEditSLA./100)/sum(n1)*100;
editSLAstd = sum(n1.*data.percEditSLAstd./100)/sum(n1)*100;
editIB = sum(n1.*data.percEditIB./100)/sum(n1)*100;
BiasedOrbit = sum(n1.*data.percBiasedOrbit./100)/sum(n1)*100;
editWet = sum(n1.*data.percEditWet./100)/sum(n1)*100;
editDry = sum(n1.*data.percEditDry./100)/sum(n1)*100;
editIono = sum(n1.*data.percEditIono./100)/sum(n1)*100;
editSsb = sum(n1.*data.percEditSsb./100)/sum(n1)*100;
editSigma0SLA = sum(n1.*data.percEditSigma0SLA./100)/sum(n1)*100;
editSigma0stdSLA = sum(n1.*data.percEditSigma0stdSLA./100)/sum(n1)*100;
editSigma0 = sum(n3.*data.percEditSigma0./100)/sum(n3)*100;
editSigma0std = sum(n3.*data.percEditSigma0std./100)/sum(n3)*100;
editSWH = sum(n2.*data.percEditSWH./100)/sum(n2)*100;
editSWHstd = sum(n2.*data.percEditSWHstd./100)/sum(n2)*100;
editWsp = sum(n4.*data.percEditWsp./100)/sum(n4)*100;
% editAll refers to the percentage of flag-valid records (over ocean and
% lakes and no Polar regions) that have been rejected by NOC criteria
editAllSLA = 100-sum(data.nrecNOC_ssha)/sum(n1)*100;
editAllSWH = 100-sum(data.nrecNOC_swh)/sum(n2)*100;
editAllSigma0 = 100-sum(data.nrecNOC_sigma0)/sum(n3)*100;
% -------------------------------------------------------------------------


% save percentage of editted data
fprintf(fid,'%.1f\n',editFlag.ssha);
fprintf(fid,'%.1f\n',editFlag.swh);
fprintf(fid,'%.1f\n',editFlag.sigma0);
fprintf(fid,'%.1f\n',editFlag.wsp);
fprintf(fid,'%.1f\n',editFlag.misp);
fprintf(fid,'%.1f\n',editSLA);
fprintf(fid,'%.1f\n',editSLAstd);
fprintf(fid,'%.1f\n',editIB);
fprintf(fid,'%.1f\n',BiasedOrbit);
fprintf(fid,'%.1f\n',editWet);
fprintf(fid,'%.1f\n',editDry);
fprintf(fid,'%.1f\n',editIono);
fprintf(fid,'%.1f\n',editSsb);
fprintf(fid,'%.1f\n',editSigma0SLA);
fprintf(fid,'%.1f\n',editSigma0stdSLA);
fprintf(fid,'%.1f\n',editSigma0);
fprintf(fid,'%.1f\n',editSigma0std);
fprintf(fid,'%.1f\n',editSWH);
fprintf(fid,'%.1f\n',editSWHstd);
fprintf(fid,'%.1f\n',editWsp);
fprintf(fid,'%.1f\n',editAllSLA);
fprintf(fid,'%.1f\n',editAllSWH);
fprintf(fid,'%.1f\n',editAllSigma0);

% set warning for corrections
if corrections.ionoNumValid(end) == 0
    fprintf(fid2,'%s\n','iono_missing');
end
if corrections.dryNumValid(end) == 0
    fprintf(fid2,'%s\n','dry_missing');
end
if corrections.wetNumValid(end) == 0
    fprintf(fid2,'%s\n','wet_missing');
end
if corrections.atmNumValid(end) == 0
    fprintf(fid2,'%s\n','atm_missing');
end
if corrections.ssbNumValid(end) == 0
    fprintf(fid2,'%s\n','ssb_missing');
end

% ----------------------------
fclose(fid);
fclose(fid2);
% ------------------------------------------------- end diagnostics
% got_to_this_point = 2222222 ;
% save /noc/users/cryo/matlab_debug.mat  data_type report_date got_to_this_point j


