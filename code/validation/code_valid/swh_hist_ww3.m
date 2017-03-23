function swh_hist_ww3(month)
% WW3 provides hourly data
%
% input : month: string in format 'mmm-yyyy' ('Mar-2015') representing the
% year up to which the validation is performed
path2 = '/noc/users/cryo/QCV_Cryo2/code/mode_mask/';
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures
% UPDATED JAN 2017
% path_ww3 = '/scratch/general/cryosat/validation_data_2/WWIII/';
% path_C2Data = '/scratch/general/cryosat/validation_data/';
path_ww3 = '/noc/mpoc/cryo/cryosat/validation_data_2/WWIII/';
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';

% ------------------------- get data from SOEST ---------------------------
%   (http://oos.soest.hawaii.edu/erddap/griddap/NWW3_Global_Best.html)

d0 = datenum(month,'mmm-yyyy');
day_date = datestr(d0+14,'yyyymmdd');
[Y,M] = datevec(d0);
ndays = eomday(Y,M);
df = d0 + ndays - 1;
Y = day_date(1:4);
M = day_date(5:6);
z = load([path_ww3 'swh_nww3_' Y M '.mat']);
fldn = fieldnames(z);
z = z.(fldn{1});
lat = double(z.latitude);
lon = double(z.longitude);
%t = double(z.time);
swh = double(z.Thgt);
%t = t/(24*3600) + datenum(1970,1,1); % seconds since 1970-01-01 (we use the whole month)
swh = squeeze(swh); %units = m
clear z
% -------------------------------------------------------------------------

d = datestr(d0:df,'yyyymmdd');
swh_c2 = [];
for i=1:ndays
    
    tmp = load([path_C2Data 'dataVal_SIR_GOP_L2_' d(i,:) '.mat']);
    swh_c2 = [swh_c2 tmp.swh]; %#ok
    
end

% ----- exclude polar regions in WWIII based on the mode mask polygons ----
load([path2 'seasMask_polar.mat']); % arctic mask

% choose seasonal mode mask
seasMonth = day_date(5:6);
seasDay = day_date(end-1:end);
seas = 2*(str2num(seasMonth)-1); %#ok
if(str2double(seasDay) < 15)
    seas = seas+1;
else
    seas = seas+2;
end
polar = polar{seas}; %#ok polar is loaded from file
lon(lon>180) = lon(lon>180) - 360;
[lon,ix] = sort(lon);
[lon,lat] = meshgrid(lon,lat);
knoPolar = isLRM(lon(:),lat(:),polar);
swh = swh(:,:,ix);
swh = reshape(swh,size(swh,1),size(swh,2)*size(swh,3));
swh(:,~knoPolar)=NaN;
% -------------------------------------------------------------------------

swh(swh>12) = NaN;
swh_c2(swh_c2>12) = NaN;
swh = swh(:);
xvalues = 0:0.25:12;
rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
mu1 = nanmean(swh);
mu2 = nanmean(swh_c2);

[h1,x1] = hist(swh,xvalues);
[h2,x2] = hist(swh_c2,xvalues);

pPos = [15 13];
pos = [12 10.5];
figure(1)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[2.2 1.7 pos])
hold on
%hb1 = bar(x1,h1,'FaceColor',rgb2,'EdgeColor','none');
hb2 = bar(x2,h2/trapz(x2,h2),0.5,'FaceColor',rgb2,'EdgeColor','none');
hb1 = plot(x1,h1/trapz(x1,h1),'LineWidth',2.5,'color',rgb1);
box on
hl = legend([hb1 hb2],{['WWIII (\mu = ' num2str(mu1,'%4.2f') ' m)'], ...
    ['CryoSat-2 (\mu = ' num2str(mu2,'%4.2f') ' m)']},'FontSize',16);
set(hl,'box','off')
set(gca,'FontSize',16)
xlim([0 12])
xlabel('SWH (m)')
ylabel('Pdf')

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_hist_swh_ww3'],'-painters','-r300') % save figure







