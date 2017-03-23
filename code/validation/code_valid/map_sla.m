function map_sla(month)
% UPDATED JAN 2017
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path4 = '/noc/mpoc/cryo/cryosat/validation_data_2/';


% path_C2Data = '/scratch/general/cryosat/validation_data/';
% path4 = '/scratch/general/cryosat/validation_data_2/';
load([path4 'gshhs_i.mat'])

% ---------------------- options ------------------------------------------
data_type = 'SIR_GOP_L2';
interp_type = 'natural'; % using boxes ('box') or natural interpolation ('natural')
dx = 1.0; % 
D0 = '20140411';
Df = datenum(month,'mmm-yyyy');
[Y,M,~] = datevec(Df);
Df = Df + eomday(Y,M) - 1;
dt = 5;
% -------------------------------------------------------------------------

% ------------------------------- read data -------------------------------
d0 = datenum(D0,'yyyymmdd');
df = Df;
d = datestr(d0:df,'yyyymmdd');
ndays = size(d,1);
% we spread the remainder of the division ndays/dt over the whole period
nweeks = fix(ndays/dt);
Rem = rem(ndays,dt);
DT = dt.*ones(nweeks,1);
DT = DT + [ones(Rem,1); zeros(nweeks-Rem,1)];
% ----------------------------------------------------------------------

lon = -200:1:200;
lat = -90:1:90;
[lat,lon] = meshgrid(lat,lon);
files = [repmat(['dataVal_' data_type '_'],ndays,1) d,...
         repmat('.mat',ndays,1)];

x = cell(nweeks,1);
y = cell(nweeks,1);
v = cell(nweeks,1);
t = cell(nweeks,1);
o = cell(nweeks,1);
for i=1:nweeks
    disp(i)
    if i ~= 1
        j0 = 1 + DT(i-1)*(i-1);
    else
        j0 = 1;
    end
    jf = j0 + DT(i) - 1;
    for j=j0:jf
        dataVal = load([path_C2Data files(j,:)]);
        xtemp = dataVal.x_ssha;
        ytemp = dataVal.y_ssha;
        vtemp = dataVal.ssha; %tides out, atm out
        ttemp = dataVal.t_ssha;
        otemp = dataVal.orbit_ssha;
        
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
    x{i} = [xr x{i} xl];
    y{i} = [y{i}(kr) y{i} y{i}(kl)];
    v{i} = [v{i}(kr) v{i} v{i}(kl)];
end
% -------------------------------------------------------------------------


if strcmp(interp_type,'box')
    z = mapbox(x,y,v,lon(:),lat(:),dx);
    z(z == 0) = NaN;
    z = reshape(z,nweeks,size(lon,1),size(lon,2));
else
    z = NaN(nweeks,length(lat(:)));
    for i=1:nweeks
        F = TriScatteredInterp(x{i}',y{i}',v{i}','natural');
        %F = scatteredInterpolant(x{i}',y{i}',v{i}','natural');
        z(i,:) = F(lon(:),lat(:));
    end
end

% restrict coordinates to (-180,180)
k = find(lon(:,1)>=-179 & lon(:,1)<=180);
z = reshape(z,size(z,1),size(lon,1),size(lon,2));
z = z(:,k,:);
lon = lon(k,:);
lat = lat(k,:);
kland = island(lon,lat,gshhs_i);
z(:,kland) = NaN; %#ok
save([path4 'ssha_c2_20140413_5days.mat'],'z','-mat')

        


