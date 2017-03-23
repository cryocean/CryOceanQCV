% validtion with HF radar (hourly data)
function HF_val(month)

%addpath(genpath('/Volumes/noc/users/cryo/QCV_Cryo2/code/'))
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures

% UPDATED JAN 2017
% path3 = '/scratch/general/cryosat/validation_data_2/HF_data/';
% path_C2Data = '/scratch/general/cryosat/validation_data/';
% path2Data = '/scratch/general/cryosat/validation_data_2/'; %gshhs
path3 = '/noc/mpoc/cryo/cryosat/validation_data_2/HF_data/';
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path2Data = '/noc/mpoc/cryo/cryosat/validation_data_2/'; %gshhs

% ---------------------- options ------------------------------------------
daten = datenum(month,'mmm-yyyy');
[Y,M] = datevec(daten);

% % % % % % % % % % % % % % % % % % % %
% ********** NB FOR MATLAB 2015 ************************
% the next line will fail if not edited in 2015 version
% ******************************************************
% % % % % % % % % % % % % % % % % % % %


Df = datestr([Y,M,eomday(Y,M),0,0,0],'yyyymmdd'); 
D0 = '20140410';

Nw = 31; %number of points to use in the along-track smoothing (odd number)
distMax = 5; % maximum distance between the satellite and the HF measurements
maxRadius = 25; % maximum radius of the area over which averaging is computed
data_set = dir([path3 'IMOS*']);
data_set = struct2cell(data_set);
data_set = data_set(1,:);
% -------------------------------------------------------------------------

XHF = cell(length(data_set),1);
YHF = cell(length(data_set),1);
UHF = cell(length(data_set),1);
VHF = cell(length(data_set),1);
THF = cell(length(data_set),1);
XO = cell(length(data_set),1);
YO = cell(length(data_set),1);
UGHF = cell(length(data_set),1);
UGC2 = cell(length(data_set),1);

d0 = datenum(D0,'yyyymmdd');
df = datenum(Df,'yyyymmdd');


for fln=1:length(data_set) % for each radar site
    disp(fln)
    disp(data_set{1})
    %------------------------- path to HF radar data --------------------------
    path2 = [path3,data_set{fln} '/'];
    fn = dir([path2 '*nc']);
    fn = struct2cell(fn);
    fn = fn(1,:)';  % hourly datafiles for radar site(fln)
    
    
    % restrict to the reference period
    %f1 = fn{1};
    %s1 = [f1(1:13) '2014'];
    
    %K=cellfun(@(x) strncmp(x,s1,17),fn);
    K=cellfun(@(x) str2double(x(14:17))>=2014,fn);
    fn = fn(K);
    
    if(fln == 4)
        load([path3 'lon_HF.mat'])
        load([path3 'lat_HF.mat'])
    else
        % read coordinates
        ncid = netcdf.open([path2 fn{1}],'NC_NOWRITE');
        
        lat = netcdf.inqVarID(ncid,'LATITUDE');
        lat = netcdf.getVar(ncid,lat,'double');
        
        lon = netcdf.inqVarID(ncid,'LONGITUDE');
        lon = netcdf.getVar(ncid,lon,'double');
        if size(lon,1) == 1 || size(lon,2) == 1
            [lat,lon] = meshgrid(lat,lon);
        end
        
        netcdf.close(ncid);
    end
    [nx,ny] = size(lon);
    XHF{fln} = lon;
    YHF{fln} = lat;
    
    
    %read HF velocity
    nhours = (df - d0 + 4)*24; % number of hours in the period of interest
    
    Uhf = zeros(nhours,nx,ny);
    Vhf = zeros(nhours,nx,ny);
    Thf = zeros(nhours,1);
    ct = 0;
    for j=1:length(fn)
        if(mod(j,1000) == 0)
            disp(j)
        end
        ncid = netcdf.open([path2 fn{j}],'NC_NOWRITE');
        
        t = netcdf.inqVarID(ncid,'TIME');
        t = netcdf.getVar(ncid,t,'double') + datenum(1950,1,1);
        if(t > d0-2 && t < df+2)
            
            u = netcdf.inqVarID(ncid,'UCUR');
            fillVal = netcdf.getAtt(ncid,u,'_FillValue','double');
            u = netcdf.getVar(ncid,u,'double');
            u(u == fillVal) = NaN;
            
            v = netcdf.inqVarID(ncid,'VCUR');
            fillVal = netcdf.getAtt(ncid,v,'_FillValue','double');
            v = netcdf.getVar(ncid,v,'double');
            v(v == fillVal) = NaN;
            
            yi = netcdf.inqVarID(ncid,'LATITUDE');
            yi = netcdf.getVar(ncid,yi,'double');
            
            xi = netcdf.inqVarID(ncid,'LONGITUDE');
            xi = netcdf.getVar(ncid,xi,'double');
            if size(xi,1) == 1 || size(xi,2) == 1
                [yi,xi] = meshgrid(yi,xi);
            end
            if ~all(size(xi) == size(lon)) % for some stations the grid changes over time
                
                Fu = TriScatteredInterp(xi(:),yi(:),u(:));
                Fv = TriScatteredInterp(xi(:),yi(:),v(:));
                u = Fu(lon,lat);
                v = Fv(lon,lat);
                
            end
            
            
            
            %         % perform linear interpolation
            %         kn = u==u;
            %         if sum(kn(:)) > 20
            %             F = TriScatteredInterp(lon(kn),lat(kn),u(kn));
            %             u = F(lon,lat);
            %
            %             kn = v==v;
            %             F = TriScatteredInterp(lon(kn),lat(kn),v(kn));
            %             v = F(lon,lat);
            %             Uhf(j,:,:) = u;
            %             Vhf(j,:,:) = v;
            %         end
            
            ct = ct+1;
            Uhf(ct,:,:) = u;
            Vhf(ct,:,:) = v;
            Thf(ct) = t;
        end
        netcdf.close(ncid); % close netcdf file

    end
    kin = Thf ~= 0;
    Uhf = Uhf(kin,:,:);
    Vhf = Vhf(kin,:,:);
    Thf = Thf(kin);
    load([path3 'U_mean_2011_2014_' num2str(fln) '.mat'])
    load([path3 'V_mean_2011_2014_' num2str(fln) '.mat'])
    Uhf = Uhf - permute(repmat(U,[1 1 sum(kin)]),[3 1 2]);
    Vhf = Vhf - permute(repmat(V,[1 1 sum(kin)]),[3 1 2]);
    UHF{fln} = U;
    VHF{fln} = V;
    kin = find(kin);
    
    disp('end reading HF radar data')
    %--------------------------------------------------------------------------
    
    
    
    
    % ------------------------------- read C2 data ----------------------------
    nw = (Nw-1)/2;
    d = datestr(d0:df,'yyyymmdd');
    %dates = d0:df;
    ndays = size(d,1);
    
    crossData = repmat(struct('x',[],'y',[],'t',[],'ssha',[],'ij',[],...
        'ir',[]),10000,1);
    xo = cell(10000,1);
    yo = cell(10000,1);
    kpass = zeros(10000,1);
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
                kd = do < distMax;
                if (sum(kd) ~= 0)
                    kmin = find(do == min(do),1);
                    dij = Distance(lon(jj,:)',lat(jj,:)',xtemp(ko(kmin)),...
                        ytemp(ko(kmin)));
                    if kmin-nw <= 0 || kmin+nw > length(ko)
                        continue
                    end
                    ts = ttemp(ko(kmin-nw:kmin+nw));
                    dt = diff(ts);
                    IJ = find(dij == min(dij));
                    if all(dt*24*3600 < 2) && ~isnan(U(jj,IJ))
                        kn = kn+1;
                        xo{kn} = xtemp(ko(do<100));
                        yo{kn} = ytemp(ko(do<100));
                        kpass(kn) = j;
                        crossData(kn).x = xtemp(ko(kmin-nw:kmin+nw));
                        crossData(kn).y = ytemp(ko(kmin-nw:kmin+nw));
                        crossData(kn).t = ts;
                        crossData(kn).ssha = vtemp(ko(kmin-nw:kmin+nw));
                        crossData(kn).ij = [jj,IJ];
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
    kempty = cellfun(@isempty,xo);
    xo(kempty) = [];
    yo(kempty) = [];
    kpass(kempty) = [];
    
    
    %read HF velocity
    ns = length(crossData);
    uhf = zeros(ns,length(kin));
    vhf = zeros(ns,length(kin));
    thf = zeros(ns,length(kin));
    for j=1:length(kin)
        for i=1:ns
            ir = crossData(i).ir;
            uhf(i,j) = nanmean(Uhf(j,ir),2);
            vhf(i,j) = nanmean(Vhf(j,ir),2);
            thf(i,j) = Thf(j);
        end
        
    end
    knan = isnan(uhf);
    knan = sum(knan,2);
    knan = find(knan ~= ns);
    uhf = uhf(knan,:);
    vhf = vhf(knan,:);
    thf = thf(knan,:);
    thf2 = thf;
    thf2(isnan(uhf)) = NaN;
    crossData = crossData(knan);
    xo = xo(knan);
    yo = yo(knan);
    kpass = kpass(knan);
    xo = cell2mat(xo');
    yo = cell2mat(yo');
    [xo,ix] = unique(xo);
    yo = yo(ix);
    THF{fln} = thf;
    XO{fln} = xo;
    YO{fln} = yo;
    
    
    
    % compute HF velocity perpendicular to altimeter track
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
        dt = abs(ti(nw+1)-thf2(i,:));
        if min(dt) < 1
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
    ut(knan) = [];
    ug(knan) = [];
    kpass(knan) = [];
    
    %add NaN to separate different passes
    kd = diff(kpass);
    kd = find(kd ~= 0);
    kd = kd + cumsum(kd == kd);
    utr = NaN(length(ut) + length(kd),1);
    ugr = NaN(length(ug) + length(kd),1);
    Kr = true(length(ug) + length(kd),1);
    Kr(kd) = false;
    utr(Kr) = ut;
    ugr(Kr) = ug;
    utr(end+1) = NaN;%#ok
    ugr(end+1) = NaN;%#ok
    
    UGHF{fln} = utr;
    UGC2{fln} = ugr;
    
end
%rho = cellfun(@(x,y) nancorrcoef(x,y),UGC2,UGHF);

% --------------------------------- plot map ------------------------------
load([path2Data 'gshhs_c.mat'])
xh = gshhs_c.lon;
yh = gshhs_c.lat;
lev = gshhs_c.level;
kl = (lev == 1 | lev == 6);
xh = xh(kl);
yh = yh(kl);
xh = cell2mat(xh);
yh = cell2mat(yh);
xh = xh(:)';
yh = yh(:)';

figure(1)
pPos = [14.2 9.4];
pos = [10.5 9.2];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo','MapLatLimit',[-40 -15],...
    'MapLonLimit',[110 155]);
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])

colormap(linspecer);
for i=1:length(UHF)
    pcolorm(YHF{i},XHF{i},UHF{i})
end

textm(-38,141,'1','FontSize',16)
textm(-36,132.7,'4','FontSize',16)
textm(-33.7,112,'2','FontSize',16)
textm(-28.7,112,'3','FontSize',16)


framem
mlabel('MLabelParallel','south','MLabelLocation',-180:10:180)
plabel('PLabelLocation',-90:5:90)
setm(ha,'FontSize',12)
tightmap
title('Regions where HF radar data are available','FontSize',12)

pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
xpos = get(gca,'position');
hc = colorbar;
set(gca,'position',xpos)
caxis([-0.4 0.4])
set(hc,'YTick',-0.4:0.2:0.4)
set(hc,'FontSize',12)
ylabel(hc,'Mean u velocity (m/s)') % add units to colorbar
cpos = get(hc,'Position');
cpos(3) = 0.6.*cpos(3);
cpos(2) = cpos(2)+0.728;
set(hc,'Position',cpos)
set(gca,'position',xpos)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_HF_map_regions'],'-r300') % save figure
close all
% -------------------------------------------------------------------------


% ----------------------- plot time series example ------------------------
u1 = cell2mat(UGHF);
u2 = cell2mat(UGC2);
region = cellfun(@length,UGHF);
y1 = u1.*100;
y2 = u2.*100;
k = find(abs(y1)>50);
y1(k) = NaN;
y2(k) = NaN;
t2 = (1:length(y1))';

pPos = [17 10];
pos = [15 8];
figure
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.7 1.5 pos])
rgb1 = [215/255 75/255 75/255];
rgb2 = [120/255 121/255 116/255];
hold on
hp1 = plot(t2,y1,'o-','color',rgb1,'LineWidth',0.8,'MarkerSize',2, ...
    'MarkerFaceColor',rgb1);
hp2 = plot(t2,y2,'o-','color',rgb2,'LineWidth',0.8,'MarkerSize',2, ...
    'MarkerFaceColor',rgb2);

hl = legend([hp1 hp2],{'HF radar','CryoSat-2'},...
    'Location','northoutside','Orientation','horizontal');
set(hl, 'Box', 'off')
set(hl, 'Color', 'none')

ylabel('Surface velocity anomaly (cm/s)','FontSize',14)
xlim([min(t2)-1 max(t2)+1])
ymin = min([y1(:);y2(:)]);
ymax = max([y1(:);y2(:)]);
ylim([ymin-10 ymax+10])

% ---------------------- ploting and naming regions -----------------------
YLIM = get(gca,'ylim');
XLIM = cumsum(region);
xLabel = {'1','2','3','4'};
xTic = [XLIM(1)/2; XLIM(1)+region(2)/2; XLIM(2)+region(3)/2; ...
    XLIM(3)+region(4)/2];

set(gca,'LineWidth',1,'FontSize',14,'XTick',xTic, ...
    'XTickLabel',xLabel,'YTick',-120:20:120)

set(gca, 'Ticklength', [0 0])
xlabel('Region number')

%region 1
x1 = (XLIM(1)+0.5);
plot([x1 ; x1], [YLIM(1) ; YLIM(2)],'color','black',...
    'lineWidth',1)

%region 2
x1 = (XLIM(2)+0.5);
plot([x1 ; x1], [YLIM(1) ; YLIM(2)],'color','black',...
    'lineWidth',1)

%region 3
x1 = (XLIM(3)+0.5);
plot([x1 ; x1], [YLIM(1) ; YLIM(2)],'color','black',...
    'lineWidth',1)

% -------------------------------------------------------------------------

box on
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out '/Fig_HF_timeSeries'],'-r300')
close all

% -------------------------------------------------------------------------





% ncid = netcdf.open('IMOS_Coffs_Harbour.nc','NC_NOWRITE');
% 
% lat = netcdf.inqVarID(ncid,'LATITUDE');
% lat = netcdf.getVar(ncid,lat,'double');
% 
% lon = netcdf.inqVarID(ncid,'LONGITUDE');
% lon = netcdf.getVar(ncid,lon,'double');
% 
% u = netcdf.inqVarID(ncid,'UCUR');
% fillVal = netcdf.getAtt(ncid,u,'_FillValue','double');
% u = netcdf.getVar(ncid,u,'double');
% u(u == fillVal) = NaN;
% 
% v = netcdf.inqVarID(ncid,'VCUR');
% fillVal = netcdf.getAtt(ncid,v,'_FillValue','double');
% v = netcdf.getVar(ncid,v,'double');
% v(v == fillVal) = NaN;
% 
% t = netcdf.inqVarID(ncid,'TIME');
% t = netcdf.getVar(ncid,t,'double');
% t = t + datenum(1950,1,1);
% 
% netcdf.close(ncid);






