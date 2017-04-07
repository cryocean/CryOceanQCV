% function read_cryo2DatavC
function read_cryo2DatavC(in_dir,out_dir,start,stop)
% Read Cryosat2 Baseline C NetCDF data and save to file as structure.
% The data are saved by day: one file per day.
%
% USAGE : read_cryo2DatavC(in_dir,out_dir,start,stop)
%
% INPUT :    in_dir : path to NetCDF data - arrange by yyyy/mm
%            out_dir : path to where data will be saved
%            start : date of the first day in format 'yyyymmdd'
%            stop : date of the last day in format 'yyyymmdd'
%
% OUTPUT : this function has no outputs
%
% author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
% change: Helen Snaith 21 Mar 2017
%         New version to read IPO/GOP baseline C NetCDF files
%
%
%% Check input arguments
% in_dir='/noc/mpoc/cryo/TDS_March_2017/SIR_IOP_L2';
% out_dir='/noc/mpoc/cryo/TDS_testout/';
% % in_dir='/noc/mpoc/cryo/TDS_March_2017/SIR_IOP_L2/';
% % out_dir='/noc/mpoc/cryo/TDS_testout/';
% % in_dir ='/Users/Shared/Data/tmp/SIR_IOP_L2/';
% % out_dir='/Users/Shared/Data/tmp/TDS_testout/';
% start = '20120603';
% stop  = '20120604';

% Input and output directories must exist and end in '/'
if ~strcmp(in_dir(end) ,'/'), in_dir  = [in_dir  '/']; end
if ~strcmp(out_dir(end),'/'), out_dir = [out_dir '/']; end
if ~exist(in_dir,'dir')
  disp(['Error: Input directory ' in_dir 'does not exist']);
  return;
end
if ~exist(out_dir,'dir')
  disp(['Error: Output directory ' out_dir 'does not exist']);
  return;
end
if ~ischar(start) || length(start)~=8 || isempty(regexp(start,'\d{8}','once'))
  disp('Error: start date must be 8 chararacter date as yyyymmdd');
  return;
end
if ~ischar(stop) || length(stop)~=8 || isempty(regexp(start,'\d{8}','once'))
  disp('Error: stop date must be 8 chararacter date as yyyymmdd');
  return;
end

%% Determine files to be read from start and stop dates
% Find timestamps from dates
t0 = datenum(start,'yyyymmdd');
tf = datenum(stop,'yyyymmdd');
if t0<datenum(2000,1,1)
  disp(['Error: start date set ' start '(' datestr(t0) ') - must be after 1/1/2000']);
  return;
end
if tf>datenum(2050,1,1)
  disp(['Error: stop date set ' stop '(' datestr(tf) ') - must be before 1/1/2050']);
  return;
end
% generate monthly directory names (yyyy/mm) from dates to read (+/1 1day)
dt = datestr((t0-1):(tf+1),'yyyy/mm/dd'); % convert date to yyyy/mm/dd string
dirs = unique(cellstr(dt(:,1:7)));        % Find unique yyyy/mm strings
dt = dt(2:end-1,:);   % remove extra +/- 1 day
dt(:,strfind(dt(1,:),'/')) = []; % change date format to yyyymmdd
dt = cellstr(dt); % convert date cell array to character array
% pre-allocate 'fln' : an [ndays x 1] cell array.
%   The j-th cell contains the paths of all files with data for the j-th day.
fln = cell(length(dt),1);

% Search all expected months for data
for i=1:length(dirs)
  if exist([in_dir dirs{i}],'dir') % skip month if there is no directory
    % Set filename search pattern based on directory name
    if strcmp(in_dir(end-10:end-1),'SIR_FDM_L2')  % L2 FDM direcotory
     not checked for FDM
        fn = dir([in_dir dirs{i} '/CS*001.DBL']); % list files in directory
      fn = struct2cell(fn);
      fn = fn(1,:);
    else % there may be several versions of the same file, we use the latest version
      fn = dir([in_dir dirs{i} '/CS*.DBL']); % directory list of /DBL files
      fn={fn(:).name};                       % save filenames to cell array
      fn2 = cellfun(@(x) x(1:end-9),fn,'unif',0); % array of filenames without _ver.DBL
      [~,ix] = unique(fn2,'last'); % find index in array of last occurence of unique names
      ix = sort(ix);
      fn = fn(ix);  % use (sorted) index to find last version of each file
    end
    % extract dates from file name - there will be 2 matches per filename
    %tf = cellfun(@(x) regexp(x,'(20\d*)T','tokens'),fn,'unif',false);
    tf = regexp(fn,'(20\d*)T','tokens'); % think this does the same thing
    tf = [tf{:}];
    tf = reshape(tf,2,[]);
    for j=1:length(fln)
      kdti = ismember([tf{1,:}],dt{j}); % does first date in filename fall on date
      kdtf = ismember([tf{2,:}],dt{j}); % does second date in filename fall on date
      kdt = or(kdti,kdtf); % if either date is true - it includes data for this day
      if sum(kdt) ~= 0
        fnj = cellfun(@(x) [in_dir dirs{i} '/' x],fn(kdt),'unif',false);
        fln{j} = [fln{j};cell2mat(fnj')]; % cell containing array of chars
      end
    end
  end
end


for i=1:length(fln) % loop over days
  if ~isempty(fln{i}) % as long as there are files for this day
    for j=1:size(fln{i},1) % loop over files
      fn = fln{i}(j,:);
      %       kf = strfind(fn,'/');
      %       f = fn(kf(end)+1:end-4);
      %       dir_in = fn(1:(kf(end)-1));
      %       fname = fn(kf(end)+1:end);
      % replace above with system independant matlab function
      [~, f, ~] = fileparts(fn);
      % read NetCDF file into structure by filename
      [v.(f),g.(f)]=rd_c2nc(fn);
      % Convert the variables read from NetCDF file to their equivalents
      v.(f) = c2_vC2vB(v.(f));
      
    end
    kp = strfind(fn,'SIR');
    prod = fn(kp(1):(kp(1)+9));
    save_cryo2DatavC(v,g,dt{i},out_dir,prod); % concatenate and save data
    clear v g % reset structures for each day
  end
end

end

function save_cryo2DatavC(v,g,day,out_dir,product)

% concatenate all variables together.
v = struct2cell(v); % cell containing structure for each file [nfiles x 1]
kempty = cellfun(@isempty,v);
v(kempty) = [];
v = cell2mat(v); % array containing structure for each file [nfiles x 1]
% Find the day for each 1Hz & each 20Hz record
kday_01 = cell(length(v),1);
kday_20 = cell(length(v),1);
n_01 = length([v(:).time]);
n_20 = length([v(:).hi_time]);
for j=1:length(v)
  tj = datenum(2000,1,1)+v(j).utc_time./(24*3600);
  kday_01{j} = ismember(cellstr(datestr(tj,'yyyymmdd')),day);
  tj = datenum(2000,1,1)+v(j).hi_utctime./(24*3600);
  kday_20{j} = ismember(cellstr(datestr(tj,'yyyymmdd')),day);
end
nrec_01 = cellfun(@sum,kday_01); % total number of 1Hz records on day
nrec_20 = cellfun(@sum,kday_20); % total number of 20Hz records on day
kday_01 = cell2mat(kday_01); % indices of 1Hz records on day
kday_20 = cell2mat(kday_20); % indices of 120Hz records on day

fldn = fieldnames(v); % Variable list
for i=1:length(fldn)  % For each variable
  v1 = [v(:).(fldn{i})]; % Generate catenated variable
  if length(v1)==n_01
    v1 = v1(kday_01);  % Subset 1Hz data to those on day
  elseif  length(v1)==n_20
    v1 = v1(kday_20);  % Subset 20Hz data to those on day
  end
  vo.(fldn{i}) = v1; %#ok
end
v0.nrec_01 = nrec_01; % numer of records per file
v0.nrec_20 = nrec_20; % numer of records per file

fn_o = [out_dir product '_' day '.mat'];
save(fn_o,'-struct','vo')
%save([out_dir 'sph_' product '_' day '.mat'],'-struct','sph')
%save([out_dir 'mph_' product '_' day '.mat'],'-struct','mph')

end
