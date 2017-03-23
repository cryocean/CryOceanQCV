not sure this is used
%colorado University provides data every 10 days.
clear all
addpath(genpath('/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/'))
path1 = '/Volumes/scratch/general/cryosat/daily_stats/';
dir_out = '/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/monthly_report/figures/';
path4 = '/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/';
path_out = '/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/GOP/';
%load([path4 'gshhs_i.mat'])

% ---------------------- options ------------------------------------------
data_type = 'SIR_GOP_L2';
vari = 'ssha'; % variable
interp_type = 'box'; % using boxes ('box') or natural interpolation ('natural')
dx = 2.0; %
dy = 2.0;
D0 = '20140501';
Df = '20150331';
% -------------------------------------------------------------------------

% ------------------------------- read data -------------------------------
prod = data_type(5:7);
d0 = datenum(D0,'yyyymmdd');
df = datenum(Df,'yyyymmdd');
d = datestr(d0:df,'yyyymmdd');
dateN = d0:df;
date0 = datenum('20140509','yyyymmdd');
dates = date0:10:df;
ndays = size(d,1);
nmonths = length(dates);

lon = -200.5:1:200.5;
lat = -90:1:90;
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
    k0 = dates(i)-5;
    kf = dates(i)+4;
    k0 = find(dateN == k0);
    kf = find(dateN == kf);
    kj = k0:kf;
    for j=1:length(kj)
        data = load([path1 files(kj(j),:)]);
        flag1 = data.(['validFlag_' vari]);
        flag2 = data.(['validNOC_' vari]);
        flag3 = data.surface_type;
        xtemp = data.lon(flag1 & flag2 & flag3 == 0);
        ytemp = data.lat(flag1 & flag2 & flag3 == 0);
        vtemp = data.(vari)(flag1 & flag2 & flag3 == 0);
        ttemp = data.tutc(flag1 & flag2 & flag3 == 0);
        otemp = data.orbits(flag1 & flag2 & flag3 == 0);
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

gmsl = zeros(nmonths,1);
lat2 = -65:1:65;
ai = cosd(lat2);
for i=1:nmonths
    Ki = true(size(lat2));
    for j=1:length(lat2)
        k = y{i}>=lat2(j)-0.5 & y{i}<lat2(j)+0.5;
        if sum(k) == 0
            Ki(j) = false;
        else
            nk = sum(k);
            wi = ai(j)/nk;
            gmsl(i) = gmsl(i) + sum(v{i}(k).*wi);
        end
    end
    gmsl(i) = gmsl(i)/sum(ai(Ki));
end
gmsl = gmsl.*100;

pause
% gmsl5 = zeros(nmonths,1);
% wi = cosd(lat2);
% for i=1:nmonths
%     Ki = true(size(lat2));
%     for j=1:length(lat2)
%         k = y{i}>=lat2(j)-0.5 & y{i}<lat2(j)+0.5;
%         if sum(k) == 0
%             Ki(j) = false;
%         else
%             gmsl5(i) = gmsl5(i) + mean(v{i}(k))*wi(j);
%         end
%     end
%     gmsl5(i) = gmsl5(i)/sum(wi(Ki));
% end
% gmsl5 = gmsl5.*100;

yy = load('/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/GOP/GMSL_CU_2014_SEAS.txt'); % GMSL Colorado University
gmsl3 = yy(:,2);
gmsl3 = gmsl3./10;
gmsl3 = gmsl3(13:end);
gmsl = gmsl-mean(gmsl(1:length(gmsl3)))+mean(gmsl3);

pause


dl = datestr(dates);
dl = cellstr(dl);
rgb1 = [215/255 75/255 75/255];
rgb2 = [120/255 121/255 116/255];
hold on
h1 = plot(gmsl3,'-o','LineWidth',2,'MarkerFaceColor',rgb1,'color',rgb1);
h2 = plot(gmsl,'-o','lineWidth',2,'color',rgb2,'MarkerFaceColor',rgb2);
hl = legend([h1,h2],'Colorado University','CryoSat-2');
set(gca,'FontSize',16)
ylabel('GMSL (cm)')
ylim([5.0 8.5])
set(gca,'XTick',2:5:length(gmsl),'XTickLabel',dl(2:5:end))
box on
xlim([0 length(gmsl)+1])
rotateXLabels(gca,45)
set(gcf,'paperPositionMode','auto')
print([path_out 'Fig_GMSL_C2_CU'],'-r300','-depsc')

    
    

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

        


