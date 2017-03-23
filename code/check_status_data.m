clear all
path_code = '/noc/users/cryo/QCV_Cryo2/code/';
addpath /noc/users/cryo/QCV_Cryo2/code/validation/code_valid/
path_return = pwd;
% UPDATED JAN 2017
root_path = '/noc/mpoc/cryo/';

% root_path = '/scratch/general/';

% Tide gauge data
[xtg,ytg,ztg,ttg] = read_tgData_ascii(335);
k = find(ztg == ztg,1,'last');
disp(['Tide gauge data: ',datestr(ttg(k),'dd-mm-yyyy')])

%buoy data
path1 = [root_path 'cryosat/validation_data_2/buoys/'];
fn = dir([path1,'*_swh']);
fn = struct2cell(fn);
fn = fn(1,:);
fn = fn(end);
yr = fn{1}(1:4);
cd([path1,fn{1}])
fn = dir;
fn = struct2cell(fn);
fn = fn(1,:);
k = cellfun(@(x) strcmp(x(1),'.'),fn);
fn(k) = [];
dt = cellfun(@(x) datenum(x,'mmm'),fn);
k = find(dt == max(dt));
if(isempty(fn))
    yr = str2num(yr)-1; %#ok
    disp(['Buoy data: ','Dec','-',num2str(yr)])
else
    disp(['Buoy data: ',fn{k},'-',yr])
end

%WWIII model
path1 = [root_path 'cryosat/validation_data_2/WWIII/'];
fn = dir([path1,'*mat']);
fn = struct2cell(fn);
fn = fn(1,:);
fn = fn{end}(end-9:end-4);
disp(['WWIII model data: ',datestr(datenum(fn,'yyyymm'),'mmm-yyyy')])

%HF radar
path3 = [root_path '/cryosat/validation_data_2/HF_data/'];

data_set = dir([path3 'IMOS*']);
data_set = struct2cell(data_set);
data_set = data_set(1,:);

path2 = [path3,data_set{1} '/'];
fn = dir([path2 '*nc']);
fn = struct2cell(fn);
fn = fn(1,:)';
fn = fn{end};
fn = fn(14:21);
disp(['HF radar data: ',datestr(datenum(fn,'yyyymmdd'),'dd-mmm-yyyy')])

%OSCAR data
path2Oscar = [root_path 'cryosat/validation_data_2/'];

%ncid = netcdf.open('oscar-third-5527cc9f9fa27.nc','NC_NOWRITE');
%ncid = netcdf.open([path2Oscar,'oscar_third.nc'],'NC_NOWRITE');
ncid = netcdf.open([path2Oscar,'oscar_vel2016.nc'],'NC_NOWRITE');

t = netcdf.inqVarID(ncid,'time');
tunits = netcdf.getAtt(ncid,t,'units');
tunits = tunits(11:end);
t = netcdf.getVar(ncid,t,'double');
t = t + datenum(tunits,'yyyy-mm-dd HH:MM:SS');
netcdf.close(ncid)
disp(['OSCAR data: ',datestr(t(end),'dd-mm-yyyy')])

%steric heights
path1 = [root_path 'cryosat/validation_data_2/EN4_TS_profiles/'];
cd(path1)
fn = dir('*.nc');
fn = struct2cell(fn);
fn = fn(1,:)';
fn = fn{end};
fn = fn(end-8:end-3);
disp(['Steric data: ',datestr(datenum(fn,'yyyymm'),'mmm-yyyy')])

%GMSL Colorado
path1 = [root_path 'cryosat/validation_data_2/'];
yy = load([path1,'GMSL_COLORADO_2014_SEAS.txt']); % GMSL Colorado University
t = yy(:,1);
t = datenum(fix(t),1,1) + rem(t,fix(t)).*365.24;
disp(['GMSL data: ',datestr(t(end),'dd-mm-yyyy')])


%GMSL Goddard
path1 = [root_path 'cryosat/validation_data_2/'];
yy = load([path1,'GMSL_GODDARD_2014_SEAS.txt']); % GMSL Goddard
t = yy(:,3);
t = datenum(fix(t),1,1) + rem(t,fix(t)).*365.24;
disp(['GMSL Goddard data: ',datestr(t(end),'dd-mm-yyyy')])


%GMSL NOAA
path1 = [root_path 'cryosat/validation_data_2/'];
yy = load([path1,'GMSL_PODAAC_2014_SEAS.txt']); % GMSL NOAA
t = yy(:,1);
t = datenum(fix(t),1,1) + rem(t,fix(t)).*365.24;
disp(['GMSL NOAA/NESDIS/STAR data: ',datestr(t(end),'dd-mm-yyyy')])

cd(path_return)









