clear all
addpath(genpath('/noc/users/fmc1q07/QCV_Cryo2/code/'))
addpath('/nerc/packages/satprogs/satmat');
satscript
data_type = 'SIR_GOP_L2';
vari = 'ssha';
dx = 0.5;
path1 = '/noc/users/fmc1q07/QCV_Cryo2/code/GOP/sla/';
d0 = '20140410';
df = '20140910';
d0 = datenum(d0,'yyyymmdd');
df = datenum(df,'yyyymmdd');
d = datestr(d0:df,'yyyymmdd');
ndays = length(d);
dt = 2;
nweeks = ndays/dt;
load([path1 'lon.txt']);
load([path1 'lat.txt']);
lon = -20:1:380;
lat = -90:1:90;
[lat,lon] = meshgrid(lat,lon);
files = [repmat(['GOP_sla_' data_type '_'],ndays,1) d repmat('.mat',ndays,1)];
data = struct([]);
x = cell(nweeks,1);
y = cell(nweeks,1);
v = cell(nweeks,1);
for i=1:nweeks
    disp(i)
    j0 = 1+dt*(i-1);
    jf = dt+dt*(i-1);
    for j=j0:jf
        data = load([path1 files(j,:)]);
        flag1 = data.(['validFlag_' vari]);
        xtemp = data.lon(flag1);
        ytemp = data.lat(flag1);
        vtemp = data.(vari)(flag1);
        flag2 = data.(['validNOC_' vari]);
        xtemp = xtemp(flag2);
        xtemp(xtemp < 0) = xtemp(xtemp < 0)+360;
        ytemp = ytemp(flag2);
        vtemp = vtemp(flag2);
        x{i} = [x{i} xtemp];
        y{i} = [y{i} ytemp];
        v{i} = [v{i} vtemp];
    end
   % kl = x{i}>0 & x{i}<=20;
   % kr = x{i}>=340 & x{i}<360;
   % xl = x{i}(kl)+360;
   % xr = x{i}(kr)-360;
   % x{i} = [xr x{i} xl];
   % y{i} = [y{i}(kr) y{i} y{i}(kl)];
   % v{i} = [v{i}(kr) v{i} v{i}(kl)];
end
disp('pause')
pause

% h = mapbox(x(1:12),y(1:12),v(1:12),lon(:),lat(:),dx);
% h(h == 0) = NaN;

z = NaN(nweeks,length(lat(:)));
for i=1:nweeks
    disp(i)
    F = scatteredInterpolant(x{i}',y{i}',v{i}','natural');
    z(i,:) = F(lon(:),lat(:));
end
k = find(lon(:,1)>=0 & lon(:,1)<=360);
z = reshape(z,size(z,1),size(lon,1),size(lon,2));
z = z(:,k,:);
lon = lon(k,:);
lat = lat(k,:);

%h1 = squeeze(nanmean(z(end-15:end,:,:)));
h1 = squeeze(nanmean(z(end-09:end,:,:)));
h1 = h1.*100;
figure
pPos = [19.5 13];
pos = [14 8];
pcolor(lon,lat,h1),shading flat
caxis([-15 15])
fillmap
colormap(linspecer(64));
box on
set(gca,'LineWidth',1.5,'FontSize',18)
ylabel('Lat')
xlabel('Lon')
h=colorbar;
set(h,'FontSize',18)
set(gca,'XTick',0:60:360)
set(gca,'YTick',-90:30:90)
xlim([0 360])
set(gca,'XTick',0:60:360)
%ylabel(h,'SWH (m)') % add units to colorbar
ylabel(h,'SLA (cm)')
%title('SLA measured by Cryosat-2, September 2014')
title('SLA measured by Cryosat-2, 20Aug-10Sep 2014')
%set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gcf,'Units', 'centimeters','OuterPosition', [30 30 pPos]);
set(gca,'Units','centimeters','Position',[2 1.5 pos])

set(gcf, 'PaperPositionMode', 'auto')
print('-depsc','SLA_Aug-Sep-2014','-painters','-r600') % save figure


