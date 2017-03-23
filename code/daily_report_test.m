probably do not use
may need path to data changing
addpath(genpath('/noc/users/cryo/QCV_Cryo2/code/'))
report_date = '20160919'; %yyyymmdd
data_type = 'SIR_FDM_L2';
in_dir = ['/scratch/general/cryosat/' data_type '/'];
path1 = '/scratch/general/cryosat/daily_stats/';
path2 = '/scratch/general/cryosat/daily_data/';
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_daily_report/figures/';

datef = datenum(report_date,'yyyymmdd');
date0 = datef -29; %30-day window
t_window = datestr(date0:datef,'yyyymmdd');

%------------------ input data for statistics -----------------------------
path4 = '/noc/users/cryo/QCV_Cryo2/code/';
load([path4 'DTU10MSS.mat']) % dtu10 MSS for the 20-Hz data
load([path4 'gshhs_i.mat']) % shoreline data

t_window2 = datestr((date0-1):(datef+1),'yyyymmdd');
fileTracks = [repmat('groundTrack_',32,1) t_window2 repmat('.mat',32,1)];
A = [];
for i=1:size(fileTracks,1)
    a = load([path4 'groundTracks/' fileTracks(i,:)]);
    A = [A a.A]; %#ok
end
% -------------------------------------------------------------------------
prod = data_type(5:7);

%% delete unused matlab structure files
dateEnd = date0-10;
dateInit = datenum(2014,04,10);
t_windowDel = datestr(dateInit:dateEnd,'yyyymmdd');
ndel = size(t_windowDel,1);
files = [repmat([data_type '_'],ndel,1) t_windowDel repmat('.mat',ndel,1)];
fileStat = [repmat(['Stats_' data_type '_'],ndel,1) t_windowDel repmat('.mat',ndel,1)];


%% Ensure data are available in matlab structure, if not process them
files = [repmat([data_type '_'],30,1) t_window repmat('.mat',30,1)];
fileStat = [repmat(['Stats_' data_type '_'],30,1) t_window repmat('.mat',30,1)];
kempty = false(1,length(files));
for i=1:size(files,1)
    if ~exist([path2 files(i,:)],'file')
        read_cryo2Data(in_dir,path2,t_window(i,:),t_window(i,:));
    end
    if i == length(files) && ~exist([path2 files(i,:)],'file')
        error('dailyReport:noData','There are no data for present day. Exiting...') % exit if no data for present day
    elseif ~exist([path2 files(i,:)],'file')
        dailyEmptyStats(t_window(i,:),data_type,path1); % create empty structure if no data
        kempty(i) = true;
    end
end

%% Read statistics in 30-day window 
data = struct([]);
corrections = struct([]);
for i=1:size(fileStat,1)
    if ~exist([path1 fileStat(i,:)],'file')
        daily_stats_v6(t_window(i,:),data_type,path2,in_dir,path1, ...
                       DTU10MSS,gshhs_i,A);
    end
    data = [data load([path1 fileStat(i,:)])]; %#ok
    tmp = load([path1 fileStat(i,:)]);
    corrections = [corrections tmp.corrections]; %#ok
end

polar = data(end).polarMask;
orbit_day = data(end).orbits;
orbB = data(end).orbitBias_ssha;
orbitUni = unique(orbit_day);
orbf = data(end).orbitLastFull;
orb0 = data(end).orbitFirstFull;
if numel(unique(orbit_day)) == 1 && orbf == 0 && orb0 == 0;
    error('dailyReport:onlyOneOrbit','There is only one imcomplete orbit. Exiting...')
end
y_swh = data(end).swh;
y_noise2DHist_inLRM = data(end).noise_ssha4scatter_inLRM;
nan_2DHist_inLRM = data(end).nan2DHist_inLRM;
if(strcmp(data_type,'SIR_IOP_L2') || strcmp(data_type,'SIR_GOP_2'))
    y_noise2DHist_outLRM = data(end).noise_ssha4scatter_outLRM;
    nan_2DHist_outLRM = data(end).nan2DHist_outLRM;
    seas = data(end).modeSeas;
    plot_mode_mask(seas,report_date) % plot seasonal mode mask
end
orb1Day = orbit_day(1);
if orbf == 0
    orbit_day(orbit_day == orbit_day(end)) = [];
end
if orb0 == 0
    orbit_day(orbit_day == orbit_day(1)) = [];
end
datai = struct2cell(data');
kempty = cellfun(@isempty,datai);
kempty = sum(~kempty,1);
kempty = (kempty == 0);
data(kempty) = [];
clear datai

perc_orbOcean = vertcat(data.perc_orbOcean);
perc_LRMOrb = vertcat(data.perc_LRMOrb);
perc_outSARinOrb = vertcat(data.perc_outSARinOrb);
perc_orbKeys = cell2mat(keys(perc_orbOcean));
perc_orbOcean = cell2mat(values(perc_orbOcean));
perc_LRMOrb = cell2mat(values(perc_LRMOrb));
perc_outSARinOrb = cell2mat(values(perc_outSARinOrb));
if any(kempty)
    orbs = perc_orbKeys(1):perc_orbKeys(end);
    otmp = ones(1,length(orbs));
    [~,io] = intersect(orbs,perc_orbKeys);
    otmp(io) = perc_orbOcean;
    perc_orbOcean = otmp;
    otmp(io) = perc_LRMOrb;
    perc_LRMOrb = otmp;
    otmp(io) = perc_outSARinOrb;
    perc_outSARinOrb = otmp;
    perc_orbKeys = orbs;
end


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
disp('done')
pause
recDate = date0:datef;
recDate(kempty) = [];
ndays = length(recDate);
ind = zeros(1,sum(data.nrec));
ind(cumsum([1 data.nrec(1:end-1)])) = 1;
ind = cumsum(ind);
recDate = recDate(ind);

nrec = data.nrec;
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
disp([length(nrecOrbit),length(perc_LRMOrb)])
perc_orbit = nrecOrbit./(nth_orbit.*perc_LRMOrb./100).*100;

