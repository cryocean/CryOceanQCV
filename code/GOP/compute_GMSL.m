not sure this is used
clear all
addpath(genpath('/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/'))
path1 = '/Volumes/scratch/general/cryosat/daily_stats/';
dir_out = '/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/monthly_report/figures/';
path4 = '/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/';
%load([path4 'gshhs_i.mat'])

% ---------------------- options ------------------------------------------
data_type = 'SIR_GOP_L2';
vari = 'ssha'; % variable
interp_type = 'box'; % using boxes ('box') or natural interpolation ('natural')
dx = 2.0; %
dy = 2.0;
D0 = '20140501';
Df = '20150130';
% -------------------------------------------------------------------------

% ------------------------------- read data -------------------------------
prod = data_type(5:7);
d0 = datenum(D0,'yyyymmdd');
df = datenum(Df,'yyyymmdd');
d = datestr(d0:df,'yyyymmdd');
ndays = size(d,1);
[Y,M] = datevec(d0:df);
Mu = unique(M);
nmonths = length(Mu);

lon = -200:2:200;
lat = -90:2:90;
[lat,lon] = meshgrid(lat,lon);
files = [repmat(['Stats_' data_type '_'],ndays,1) d repmat('.mat',ndays,1)];
data = struct([]);
x = cell(nmonths,1);
y = cell(nmonths,1);
v = cell(nmonths,1);
x2 = cell(nmonths,1);
y2 = cell(nmonths,1);
v2 = cell(nmonths,1);
t = cell(nmonths,1);
o = cell(nmonths,1);
for i=1:nmonths
    disp(i)
    kj = find(M == Mu(i));
    for j=1:length(kj)
        data = load([path1 files(kj(j),:)]);
        flag1 = data.(['validFlag_' vari]);
        flag2 = data.(['validNOC_' vari]);
        xtemp = data.lon(flag1 & flag2);
        ytemp = data.lat(flag1 & flag2);
        vtemp = data.(vari)(flag1 & flag2);
        ttemp = data.tutc(flag1 & flag2);
        otemp = data.orbits(flag1 & flag2);
        x{i} = [x{i} xtemp];
        y{i} = [y{i} ytemp];
        v{i} = [v{i} vtemp];
        t{i} = [t{i} ttemp];
        o{i} = [o{i} otemp];
    end
    kl = x{i}<=-160;
    kr = x{i}>=160;
    xl = x{i}(kl)+360;
    xr = x{i}(kr)-360;
    x2{i} = [xr x{i} xl];
    y2{i} = [y{i}(kr) y{i} y{i}(kl)];
    v2{i} = [v{i}(kr) v{i} v{i}(kl)];
end
% -------------------------------------------------------------------------
disp('done')
lat2 = -67.5:5:67.5;
gmsl = zeros(nmonths,1);
tmp = zeros(length(lat2),1);
for i=1:nmonths
    w = 0;
    for j=1:length(lat2)
        k = y{i}>=lat2(j)-2.5 & y{i}<lat2(j)+2.5;
        tmp(j) = sum(v{i}(k).*cosd(lat2(j)));
        w = w + cosd(lat2(j))*sum(k);
    end
    gmsl(i) = sum(tmp./w);
end
    
    

% interpolation
z = NaN(nmonths,length(lat(:)));
for i=1:nmonths
    %F = scatteredInterpolant(x2{i}',y2{i}',v2{i}','natural');
    F = TriScatteredInterp(x2{i}',y2{i}',v2{i}','natural');
    z(i,:) = F(lon(:),lat(:));
end

% restrict coordinates to (-180,180)
k = find(lon(:,1)>=-180 & lon(:,1)<=180);
z = reshape(z,size(z,1),size(lon,1),size(lon,2));
z = z(:,k,:);
lon = lon(k,:);
lat = lat(k,:);
z = z(:,:);
lon = lon(:)';
lat = lat(:)';
k65 = abs(lat) <= 65;
lon = lon(k65);
lat = lat(k65);
z = z(:,k65);
gmsl2 = zeros(nmonths,1);
for i=1:nmonths
    knan = z(i,:) == z(i,:);
    W = cosd(lat(knan));
    gmsl2(i) = sum(bsxfun(@times,z(i,knan),W))/sum(W);
end

        


