not sure this is used
clear all
addpath(genpath('/noc/users/fmc1q07/QCV_Cryo2/code/'))
path1 = '/scratch/general/cryosat/daily_stats/';
dir_out = '/noc/users/fmc1q07/QCV_Cryo2/code/monthly_report/figures/';
path4 = '/noc/users/fmc1q07/QCV_Cryo2/code/';
load([path4 'gshhs_i.mat'])

% ---------------------- options ------------------------------------------
data_type = 'SIR_GOP_L2';
vari = 'wsp'; % variable
interpol = 0; % tracks (0) or interpolated (1) data
interp_type = 'natural'; % using boxes ('box') or natural interpolation ('natural')
dx = 2.0; % 
D0 = '20140915';
Df = '20140917';
dt = 1;
lonLim = [-90 20]; % [0 360] will give Pacific view
latLim = [10 50];
stat = 'max'; % 'mean', 'median', or 'max'. For tracks 'max' plot the
              % highest values on top
ts = cellstr(['20140915';'20140917']); % leave empty if not used
%ts = cellstr('20140915');
% -------------------------------------------------------------------------

% ------------------------------- read data -------------------------------
prod = data_type(5:7);
d0 = datenum(D0,'yyyymmdd');
df = datenum(Df,'yyyymmdd');
d = datestr(d0:df,'yyyymmdd');
ndays = size(d,1);
nweeks = ndays/dt;
lon = -200.5:1:200.5;
lat = -90:1:90;
[lat,lon] = meshgrid(lat,lon);
files = [repmat(['Stats_' data_type '_'],ndays,1) d repmat('.mat',ndays,1)];
data = struct([]);
x = cell(nweeks,1);
y = cell(nweeks,1);
v = cell(nweeks,1);
t = cell(nweeks,1);
o = cell(nweeks,1);
for i=1:nweeks
    disp(i)
    j0 = 1+dt*(i-1);
    jf = dt+dt*(i-1);
    for j=j0:jf
        data = load([path1 files(j,:)]);
        flag1 = data.(['validFlag_' vari]);
        flag2 = data.(['validNOCPolarIncluded_' vari]);
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
    if interpol == 1
        kl = x{i}<=-160;
        kr = x{i}>=160;
        xl = x{i}(kl)+360;
        xr = x{i}(kr)-360;
        x{i} = [xr x{i} xl];
        y{i} = [y{i}(kr) y{i} y{i}(kl)];
        v{i} = [v{i}(kr) v{i} v{i}(kl)];
    end
end
% -------------------------------------------------------------------------

% ------------------------- plot parameters -------------------------------
if strcmp(vari,'ssha')
    j = 1;
elseif strcmp(vari,'swh')
    j=2;
elseif strcmp(vari,'wsp')
    j=4;
end
units = [{'SSH anomaly (cm) '};{'SWH (m) '};{'Sigma0 (dB) '};...
    {'Altimeter wind speed (m/s)'}; ...
    {'Square of off-nadir angle (10^{-2} deg^2) '}];
titles = {'SSH anomaly','SWH','Sigma0','Altimeter wind speed', ...
    'Square of off-nadir angle'};
cytick = [5,1,1,5,1];
scale = [100,1,1,1,100];
cmax = [20,10,50,35,5];
cmin = [-20,0,0,0,1];
titl = [titles{j} ' measured by CryoSat-2, September 2014'];
pPos = [13.9 9.2];
pos = [10.5 9];
% -------------------------------------------------------------------------

if interpol == 0 % no interpolation
    x = [x{:}];
    y = [y{:}];
    v = [v{:}].*scale(j);
    t = [t{:}];
    o = [o{:}];
    if strcmp(stat,'max')
        [v,iv] = sort(v,2,'descend');
        y = y(iv);
        x = x(iv);
        t = t(iv);
        o = o(iv);
    end
    if ~isempty(ts)
        kdays = ismember(cellstr(datestr(t,'yyyymmdd')),ts);
        x = x(kdays);
        y = y(kdays);
        v = v(kdays);
        o = o(kdays);
        t = t(kdays);
    end
    plot_map(1,x,y,v,cmax(j),cmin(j),cytick(j),units{j},titl,pPos,pos, ...
        interpol,lonLim,latLim);
    
else % with interpolation
    if interpol == 1 && strcmp(interp_type,'box')
        z = mapbox(x,y,v,lon(:),lat(:),dx);
        z(z == 0) = NaN;
        z = reshape(z,nweeks,size(lon,1),size(lon,2));
    elseif interpol == 1
        z = NaN(nweeks,length(lat(:)));
        for i=1:nweeks
            F = scatteredInterpolant(x{i}',y{i}',v{i}','natural');
            z(i,:) = F(lon(:),lat(:));
        end
    end
    
    % restrict coordinates to (-180,180)
    k = find(lon(:,1)>=-180 & lon(:,1)<=180);
    z = reshape(z,size(z,1),size(lon,1),size(lon,2));
    z = z(:,k,:);
    lon = lon(k,:);
    lat = lat(k,:);
    kland = island(lon,lat,gshhs_i);
    if strcmpi(viewm,'Pacific')
        lon(lon<0) = lon(lon<0)+360;
        [~,IX] = sort(lon(:,1));
        lon = lon(IX,:);
        lat = lat(IX,:);
        z = z(:,IX,:);
        kland = reshape(kland,size(lon,1),size(lon,2));
        kland = kland(IX,:);
        kland = kland(:);
    end
    
    if strcmp(stat,'mean')
        z = squeeze(nanmean(z)).*scale(j);
    elseif strcmp(stat,'median')
        z = squeeze(nanmedian(z)).*scale(j);
    else
        z = squeeze(nanmax(z)).*scale(j);
    end
    z(kland) = NaN;
   
    plot_map(1,lon,lat,z,cmax(j),cmin(j),cytick(j),units{j},titl,pPos,pos, ...
        interpol,lonLim,latLim);
end

if interpol == 0
    intpol = 'no';
else
    intpol = 'yes';
end

% ------------------------hurricane track ---------------------------------
xy = load('/noc/users/fmc1q07/QCV_Cryo2/code/GOP/edouard_hurricane.txt');
plotm(xy(:,1),xy(:,2),'color',[0 0 0])
times = {'15/21','17/03','Hurricane Edouard'};
textm([28.0 33.50 xy(1,1)],[-56.5 -56.4 xy(1,2)],times)
% -------------------------------------------------------------------------
set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out prod '_Fig_' vari '_interp_' intpol, ...
    '_' stat '_' D0 '_' Df],'-painters','-r300') % save figure
close all

% -------------------- plot along track time series -----------------------
if interpol == 0 % no interpolation
    cytick = [5,2,1,5,1];
    cmax = [20,12,50,35,5];
    cmin = [-20,0,0,0,1];
    khurr = x>-60 & x<-40 & y>20 & y<45;
    x = x(khurr);
    y = y(khurr);
    v = v(khurr);
    o = o(khurr);
    t = t(khurr);
    [ou,~] = unique(o);
    for i=1:length(ou)
        ko = o == ou(i);
        xtmp = x(ko);
        ytmp = y(ko);
        vtmp = v(ko);
        [ytmp,ix] = sort(ytmp);
        vtmp = vtmp(ix);
        
        pPos = [10 10];
        pos = [10 10];
        figure(1)
        set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
        set(gca,'Units','centimeters','Position',[1.1 0 pos], ...
            'outerPosition',[1.1 0 10 10])
        plot(vtmp,ytmp,'LineWidth',2)
        xlim([cmin(j) cmax(j)])
        set(gca,'FontSize',12,'XTick',cmin(j):cytick(j):cmax(j),'YTick',20:5:45)
        yt = get(gca,'Ytick');
        ytl = cellstr(get(gca,'YTickLabel'));
        xl = get(gca,'Xlim');
        set(gca,'YTickLabel',[])
        w = xl(2)-xl(1);
        offs = 2.45*w/12;
        for k=1:numel(yt);
            text(xl(1)-offs,yt(k)+0.2,[ytl{k},'^{\circ}','N'],'FontSize',12);
        end
        xlabel(units{j})
        box on
        set(gcf, 'PaperPositionMode', 'manual')
        print('-depsc',[dir_out prod '_Fig_' vari '_orbit_' num2str(ou(i))], ...
            '-painters','-r300') % save figure
        close all
    end
end
        


