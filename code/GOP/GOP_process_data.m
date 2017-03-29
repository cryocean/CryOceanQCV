if this is used needs checking as although GOP loads FDM data
function GOP_process_data() 
addpath(genpath('/noc/users/cryo/QCV_Cryo2/code/'))
%report_date = '20140913'; %yyyymmdd
data_type = 'SIR_FDM_L2';
in_dir = ['/scratch/general/cryosat/' data_type '/'];
path1 = '/scratch/general/cryosat/daily_stats/';
path2 = '/scratch/general/cryosat/daily_data/';

date0 = datenum('20140410','yyyymmdd'); %GOP data starts on that day
%date0 = datenum('20140801','yyyymmdd');
%datef = datenum(date,'dd-mmm-yyyy')-35; % GOP is released with a latency of ~30days
datef = datenum('20140509','yyyymmdd');
t_window = datestr(date0:datef,'yyyymmdd');
ndays = size(t_window,1); % number of days in month

%------------------ input data for statistics -----------------------------
path4 = '/noc/users/cryo/QCV_Cryo2/code/';
load([path4 'DTU10MSS.mat']) % dtu10 MSS for the 20-Hz data
load([path4 'gshhs_h.mat']) % shoreline data

t_window2 = datestr((date0-1):(datef+1),'yyyymmdd');
fileTracks = [repmat('groundTrack_',ndays+2,1) t_window2 repmat('.mat',ndays+2,1)];
A = [];
for i=1:size(fileTracks,1)
    a = load([path4 'groundTracks/' fileTracks(i,:)]);
    A = [A a.A]; %#ok
end
% -------------------------------------------------------------------------

%% Store data in Matlab Structure
nfiles = size(t_window,1);
files = [repmat([data_type '_'],nfiles,1) t_window repmat('.mat',nfiles,1)];
kempty = false(1,size(files,1));
for i=1:size(files,1)
    if ~exist([path2 files(i,:)],'file')
        read_cryo2Data(in_dir,path2,t_window(i,:),t_window(i,:));
    end
    if i == size(files,1) && ~exist([path2 files(i,:)],'file')
        error('dailyReport:noData','There are no data for present day. Exiting...') % exit if no data for present day
    elseif ~exist([path2 files(i,:)],'file')
        dailyEmptyStats(t_window(i,:),data_type,path1); % create empty structure if no data
        kempty(i) = true;
    end
end

%% Compute statistics
files = [repmat(['Stats_' data_type '_'],nfiles,1) t_window repmat('.mat',nfiles,1)];
data = struct([]);
for i=1:size(files,1) 
    if ~exist([path1 files(i,:)],'file')
        daily_stats_v6(t_window(i,:),data_type,path2,in_dir,path1, ...
                       DTU10MSS,gshhs_h,A);
    end
    data = [data load([path1 files(i,:)])]; %#ok
end
