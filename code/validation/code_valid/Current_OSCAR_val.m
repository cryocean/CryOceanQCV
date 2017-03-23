% validtion with OSCAR radar
function Current_OSCAR_val(month,region)
% input : month: string in format 'mmm-yyyy' ('Mar-2015')
%         region : string denoing the region ('Atlantic' or 'Pacific') 

%addpath(genpath('/Volumes/noc/users/cryo/QCV_Cryo2/code/'))
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/';
% UPDATE HERE JAN 2017
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path2Oscar = '/noc/mpoc/cryo/cryosat/validation_data_2/';
% path_C2Data = '/scratch/general/cryosat/validation_data/';
% path2Oscar = '/scratch/general/cryosat/validation_data_2/';


% ---------------------- options ------------------------------------------
nw = 25; %number of points to use in the along-track smoothing (odd number)
distMax = 10; % maximum distance between the satellite and the HF measurements
maxRadius = 10; % maximum radius of the area over which averaging is computed

d0 = datenum(month,'mmm-yyyy');
[Y,M,~]  = datevec(d0);
df = d0 + eomday(Y,M) - 1;
% -------------------------------------------------------------------------



%ncid = netcdf.open('oscar-third-5527cc9f9fa27.nc','NC_NOWRITE');
oscar_fn = ['oscar_vel' num2str(Y) '.nc'];

ncid = netcdf.open([path2Oscar,oscar_fn],'NC_NOWRITE');




lat = netcdf.inqVarID(ncid,'latitude');
lat = netcdf.getVar(ncid,lat,'double');

lon = netcdf.inqVarID(ncid,'longitude');
lon = netcdf.getVar(ncid,lon,'double');

t = netcdf.inqVarID(ncid,'time');
tunits = netcdf.getAtt(ncid,t,'units');
tunits = tunits(11:end);
t = netcdf.getVar(ncid,t,'double');
t = t + datenum(tunits,'yyyy-mm-dd HH:MM:SS');

% u = netcdf.inqVarID(ncid,'u_anom');
% u = squeeze(netcdf.getVar(ncid,u,'double'));
u = netcdf.inqVarID(ncid,'u');
u = squeeze(netcdf.getVar(ncid,u,'double'));

v = netcdf.inqVarID(ncid,'v');
v = squeeze(netcdf.getVar(ncid,v,'double'));

netcdf.close(ncid)

kt = find(t>=d0 & t<=df);
u = u(:,:,kt);
v = v(:,:,kt);
t = t(kt);

%select a region
if strcmp(region,'Atlantic')
    %Atlantic box
    y0 = 20;
    yf = 40;
    kx = lon > 315 & lon < 325;
    ky = lat > y0 & lat < yf;
    
elseif strcmp(region,'Pacific')
    
    %Pacific box
    y0 = 20;
    yf = 40;
    kx = lon > 220 & lon < 230;
    ky = lat > y0 & lat < yf;
end

lon = lon(kx);
lat = lat(ky);
u = u(kx,ky,:);
v = v(kx,ky,:);

lon(lon>180) = lon(lon>180) - 360;
[lat,lon] = meshgrid(lat,lon);
[nx,ny] = size(lon);

u = reshape(u,nx*ny,size(u,3));
v = reshape(v,nx*ny,size(v,3));


% ------------------------------- read C2 data ----------------------------
nw = (nw-1)/2;
d = datestr(d0:df,'yyyymmdd');
ndays = size(d,1);

crossData = repmat(struct('x',[],'y',[],'t',[],'ssha',[],'ij',[],...
    'ir',[]),10000,1);
kn = 0;
for i=1:ndays
    
    dataVal = load([path_C2Data 'dataVal_SIR_GOP_L2_' d(i,:) '.mat']);
    xtemp = dataVal.x_ssha;
    ytemp = dataVal.y_ssha;
    vtemp = dataVal.ssha + dataVal.atm; %tides out, atm in
    ttemp = dataVal.t_ssha;
    otemp = dataVal.orbit_ssha;
    
    
    ou = unique(otemp);
    for j=1:length(ou)
        ko = find(otemp == ou(j));
        
        for jj=1:nx
            xj = lon(jj,:)';
            yj = lat(jj,:)';
            [yj,IX] = sort(yj);%#ok
            xj = xj(IX);
            do = distance2Coast(xtemp(ko)',ytemp(ko)',xj,yj,100)./1000;
            kd = find(do < distMax,1);
            if ~isempty(kd)
                kmin = find(do == min(do),1);
                dij = Distance(lon(jj,:)',lat(jj,:)',xtemp(ko(kmin)),...
                    ytemp(ko(kmin)));                 
                if kmin-nw <= 0 || kmin+nw > length(ko)
                    continue
                end
                ts = ttemp(ko(kmin-nw:kmin+nw));
                dt = diff(ts);
                if all(dt*24*3600 < 2)
                    kn = kn+1;
                    crossData(kn).x = xtemp(ko(kmin-nw:kmin+nw));
                    crossData(kn).y = ytemp(ko(kmin-nw:kmin+nw));
                    crossData(kn).t = ts;
                    crossData(kn).ssha = vtemp(ko(kmin-nw:kmin+nw));
                    crossData(kn).ij = [jj,find(dij == min(dij))];
                    ij = crossData(kn).ij;
                    dr = Distance(lon(:),lat(:),lon(ij(1),ij(2)),...
                            lat(ij(1),ij(2)))./1000;
                    crossData(kn).ir = find(dr < maxRadius);
                end
            end
        end
                
    end
end
kempty = arrayfun(@(X) isempty(X.ssha),crossData);
crossData = crossData(~kempty); % remove empty structures

% -------------------------------------------------------------------------



% -------------------------- prepare OSCAR data ---------------------------
ns = length(crossData);
uhf = zeros(ns,size(u,2));
vhf = zeros(ns,size(u,2));
thf = zeros(ns,size(u,2));
for j=1:size(u,2)
    for i=1:ns
        ir = crossData(i).ir;
        uhf(i,j) = nanmean(u(ir,j),1);
        vhf(i,j) = nanmean(v(ir,j),1);
        thf(i,j) = t(j);
    end
    
end

knan = isnan(uhf);
knan = sum(knan,2);
knan = find(knan ~= ns);
uhf = uhf(knan,:);
vhf = vhf(knan,:);
thf = thf(knan,:);
crossData = crossData(knan);
% -------------------------------------------------------------------------


% ----------- compute velocity perpendicular to altimeter track -----------
ns = size(uhf,1);

ut = NaN(ns,1);
ug = NaN(ns,1);
tg = NaN(ns,1);
X = NaN(ns,1);
Y = NaN(ns,1);
for i=1:ns
    xi = crossData(i).x;
    yi = crossData(i).y;
    ti = crossData(i).t;
    vi = crossData(i).ssha;
    tg(i) = ti(nw+1);
    X(i) = xi(nw+1);
    Y(i) = yi(nw+1);
    [yi,iy] = sort(yi);
    xi = xi(iy);
    ti = ti(iy);
    vi = vi(iy);
    vxy = [yi(end)-yi(1); -(xi(end)-xi(1))]; %vector normal to the track
    vxy = vxy./norm(vxy); % normalize vector for scalar product
    f = 2*7.2921e-5*sind(yi(nw+1));
    dt = abs(ti(nw+1)-thf(i,:));
    if min(dt) < 5
        uti = dot([uhf(i,:);vhf(i,:)],repmat(vxy,[1,size(uhf,2)]),1); % projection of wind
        ut(i) = interp1(thf(i,:),uti,ti(nw+1));
        
        % compute geostrophic velocity
        dx = Distance(xi(nw),yi(nw),xi(nw+1),yi(nw+1));
        [cn, ~] = genweights(nw,nw+1,0.9434);
        npq = (-nw:nw);
        npq(npq == 0) = [];
        ct = 0;
        Dh = 0;
        for im=npq
            ct = ct+1;
            Dh = Dh + cn(ct)*(vi(nw+1+im)-vi(nw+1))/(im*dx);
        end
        %dx = Distance(xi(nw),yi(nw),xi(nw+2),yi(nw+2));
        %dh = nanmean(vi(nw) - vi(nw+2));
        ug(i) = -9.81/f*Dh;
        %ug(i) = 9.81/f*dh/dx;
    end
    
end

knan=isnan(ug) | isnan(ut);
ut(knan)=NaN;
ug(knan)=NaN;


% plot as a function of latitude
[Y,ix] = sort(Y);
ut = ut(ix);
ug = ug(ix);

Y2 = ((y0+0.25):0.5:(yf+0.25))';
u1 = NaN(length(Y2),1);
u2 = NaN(length(Y2),1);
for i=1:length(Y2)
    k = find(Y > Y2(i)-0.25 & Y <= Y2(i)+0.25);
    u1(i) = nanmean(ut(k));
    u2(i) = nanmean(ug(k));
end
nancorrcoef(u1,u2)
%rms = sqrt(nanmean((u1-u2).^2));


rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
u1 = u1.*100;
u2 = u2.*100;
hold on
h1 = plot(u1,Y2,'color',rgb1,'lineWidth',2);
h2 = plot(u2,Y2,'color',rgb2,'lineWidth',2);
set(gca,'FontSize',16)
ylabel('Latitude')
xlabel('Geostrophic velocity anomaly (cm/s)')
box on
ylim([y0-1 yf+1])
hl = legend([h1,h2],'OSCAR','CryoSat-2',...
    'Location','northoutside','Orientation','horizontal');
YLIM = get(gca,'Ylim');
YLIM = YLIM(1):YLIM(2);
plot(zeros(size(YLIM)),YLIM,'--','color','black','lineWidth',1)
set(hl,'box','off')

set(gcf, 'PaperPositionMode', 'auto')
print('-depsc',[dir_out,'Fig_C2_OSCAR_',region],'-r300')
close all

