function read_cryo2DatavC
%function read_cryo2DatavC(in_dir,out_dir,start,stop)
% read binary Cryosat2 data and save them in structure. The data is saved
% by day: one file per day. 
%
% USAGE : read_cryo2Data(in_dir,out_dir,start,stop)
%
% INPUT :    in_dir : path to binary data
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
%in_dir='/Volumes/scratch/general/cryosat/SIR_FDM_L2/';
%out_dir='/Users/fmc1q07/code/read_data/daily_data/';
in_dir='/noc/mpoc/cryo/TDS_March_2017/SIR_IOP_L2/';
out_dir='/noc/mpoc/cryo/TDS_testout/';
start = '20120603';
stop  = '20120603';


%from start and stop dates determine files to be read
t0 = datenum(start,'yyyymmdd');
tf = datenum(stop,'yyyymmdd');
%generate monthly directory names from dates to read
dt = datestr((t0-1):(tf+1),'yyyy/mm/dd');
dirs = unique(cellstr(dt(:,1:7)));
dt = dt(2:end-1,:);
dt(:,strfind(dt(1,:),'/')) = [];
dt = cellstr(dt);
fln = cell(length(dt),1);
% 'fln' is an [ndays x 1] cell array. The j-th cell contains the paths of
% all files with data for the j-th day.
for i=1:length(dirs)
    if ~exist([in_dir dirs{i}],'dir')
        continue
    end
    if strcmp(in_dir(end-10:end-1),'SIR_FDM_L2')
        fn = dir([in_dir dirs{i} '/CS*001.DBL']); %list files in directory
        fn = struct2cell(fn);
        fn = fn(1,:);
    else % there may be several versions of the same file, we use the latest version
        fn = dir([in_dir dirs{i} '/CS*.DBL']); 
        fn = struct2cell(fn);
        fn = fn(1,:);
        fn2 = cellfun(@(x) x(1:end-9),fn,'unif',0);
        [~,ix] = unique(fn2,'last');
        ix = sort(ix);
        fn = fn(ix);
    end
    % extract dates from file name
    tf = cellfun(@(x) regexp(x,'(20\d*)T','tokens'),fn,'unif',false);
    tf = [tf{:}];
    tf = reshape(tf,2,length(tf)/2);
    for j=1:length(fln)
        kdti = ismember([tf{1,:}],dt{j});
        kdtf = ismember([tf{2,:}],dt{j});
        kdt = or(kdti,kdtf);
        if sum(kdt) == 0
            continue
        end
        fnj = cellfun(@(x) [in_dir dirs{i} '/' x],fn(kdt),'unif',false);
        fln{j} = [fln{j};cell2mat(fnj')]; % cell containing array of chars
    end
end


for i=1:length(fln) % loop over days
    if isempty(fln{i}) % if there are no files, continue...
        continue
    end
    for j=1:size(fln{i},1) % loop over files
        fn = fln{i}(j,:);
        kf = strfind(fn,'/');
        f = fn(kf(end)+1:end-4);
        dir_in = fn(1:(kf(end)-1));
        fname = fn(kf(end)+1:end);
% -- replace with NetCDF routines
        %[v.(f),sph.(f),mph.(f)]=rd_esa_fmc_v2(dir_in,fname,1);
        [v.(f),g.(f)]=rd_nc([dir_in '/' fname]);
    end
    kp = strfind(fn,'SIR');
    prod = fn(kp(1):(kp(1)+9));
    save_cryo2DatavC(v,g,dt{i},out_dir,prod); % concatenate and save data
    clear v sph mph
end

end


function save_cryo2DatavC(v,g,day,out_dir,product) %#ok

% Convert the variables read from NetCDF file to their equivalents

% concatenate all variables together.
v = struct2cell(v); % cell containing structure for each file [nfiles x 1]
kempty = cellfun(@isempty,v);
v(kempty) = [];
v = cell2mat(v); % array containing structure for each file [nfiles x 1]
fldn = fieldnames(v);
for i=1:length(fldn)
    v1 = [v(:).(fldn{i})];
    v2 = [v1(:).flag_mcd];
    v1 = rmfield(v1,'flag_mcd');
    names = fieldnames(v1);
    names2 = fieldnames(v2);
    kday = cell(length(v1),1);
    for j=1:length(v1)
        tj = v1(j).utc_time;
        tj = datenum(2000,1,1)+tj./(24*3600);
        kday{j} = ismember(cellstr(datestr(tj,'yyyymmdd')),day);
    end
    nrec = cellfun(@sum,kday);
    kday = cell2mat(kday);
    v1 = cellfun(@(f) {vertcat(v1.(f))'},names); %cell [nfields x 1]
    v1 = cellfun(@(f) f(:,kday),v1,'unif',false);
    v2 = cellfun(@(f) {vertcat(v2.(f))'},names2);
    v2 = cellfun(@(f) f(:,kday),v2,'unif',false);
    names2 = cellfun(@(f) ['flag_mcd_',f],names2,'unif',false); % change name of flag_mcd fields
    v1 = cell2struct([v1;v2],[names;names2]);
    v1.nrec = nrec; % numer of records per file
    vo.(fldn{i}) = v1; %#ok
end

fn_o = [out_dir product '_' day '.mat'];
save(fn_o,'-struct','vo')
save([out_dir 'sph_' product '_' day '.mat'],'-struct','sph')
save([out_dir 'mph_' product '_' day '.mat'],'-struct','mph')

end
