%colorado University provides data every 10 days.
function GMSL_computation(month)
addpath(genpath('/noc/users/cryo/QCV_Cryo2/code/'))
addpath /nerc/packages/satprogs/satmat
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures

% UPDATED JAN 2017
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path2GMSL = '/noc/mpoc/cryo/cryosat/validation_data_2/';



% path_C2Data = '/scratch/general/cryosat/validation_data/';
% path2GMSL = '/scratch/general/cryosat/validation_data_2/';
%load([path4 'gshhs_i.mat'])

% ---------------------- options ------------------------------------------
data_type = 'SIR_GOP_L2';
D0 = '20140410';
Df = datenum(month,'mmm-yyyy');
[Y,M,~] = datevec(Df);
Df = Df + eomday(Y,M) - 1;
% -------------------------------------------------------------------------

% ------------------------------- read data -------------------------------
d0 = datenum(D0,'yyyymmdd');
df = Df;
d = datestr(d0:(df+4),'yyyymmdd');
dateN = d0:(df+4);
date0 = datenum('20140415','yyyymmdd');
dates = date0:10:df;
ndays = size(d,1);
nmonths = length(dates);

% lon = -200.5:1:200.5;
% lat = -90:1:90;
%[lat,lon] = meshgrid(lat,lon);
files = [repmat(['dataVal_' data_type '_'],ndays,1) d repmat('.mat',ndays,1)];
x = cell(nmonths,1);
y = cell(nmonths,1);
v = cell(nmonths,1);
x2 = cell(nmonths,1);
y2 = cell(nmonths,1);
v2 = cell(nmonths,1);
t = cell(nmonths,1);
o = cell(nmonths,1);
for i=1:nmonths; % each month is a 10 day period
    disp(i)
    k0 = dates(i)-5; % start day
    kf = dates(i)+4; % 10 days later
    k0 = find(dateN == k0);
    kf = find(dateN == kf);
    kj = k0:kf; 
    for j=1:length(kj) % for each day
        data = load([path_C2Data files(kj(j),:)]); % load data for that day
        xtemp = data.x_ssha; % lon?
        ytemp = data.y_ssha; % lat?
        vtemp = data.ssha; % ssha?
        ttemp = data.t_ssha; % time?
        otemp = data.orbit_ssha; % orbit number?
        
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
lat2 = -65:1:65; % set latitude range
ai = cosd(lat2); % weighting for latitude (more obs at poles)
for i=1:nmonths % for each month
    Ki = true(size(lat2)); % logic for each latitude band
    for j=1:length(lat2) % for each latitude band
        k = y{i}>=lat2(j)-0.5 & y{i}<lat2(j)+0.5; % index that latitude is within the given latitude band
        if sum(k) == 0
            Ki(j) = false; % no obs in that latitude band set to false
        else % there is at least one obs in that lat band during the time period
            nk = sum(k); % number of observations
            wi = ai(j)/nk; % weighting 
            gmsl(i) = gmsl(i) + sum(v{i}(k).*wi);
        end
    end
    gmsl(i) = gmsl(i)/sum(ai(Ki));
end
gmsl = gmsl.*100;
% save /noc/users/cryo/QCV_Cryo2/code/validation/code_valid/testdata.mat gmsl y v
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

yy = load([path2GMSL,'GMSL_COLORADO_2014_SEAS.txt']); % GMSL Colorado University
%yy = load([path2GMSL,'GMSL_PODAAC_2014_SEAS.txt']); % GMSL PODAAC
gmsl3 = yy(:,2);
ty = yy(:,1);
ty = datenum(fix(ty),1,1) + rem(ty,fix(ty)).*365.24;
gmsl3 = gmsl3./10; % convert from mm to cm
gmsl3 = interp1(ty,gmsl3,dates); % interpolate UCB data to same dates as C2

% changed by CB 13 Jan 2017 as I don't undertsand how the code works, not
% to say it doesn't
% need to find where have valid data for both GMSL from C2 and UCB
valid_ucb_gmsl = vec(~isnan(gmsl3)) ;
valid_c2_gmsl = vec(~isnan(gmsl)) ;


mean_ucb = mean(gmsl3(valid_c2_gmsl & valid_ucb_gmsl)) ; % mean UCB GMSL

% change mean of GMSL from C2 to have same mean as UCB
gmsl = gmsl - mean(gmsl(valid_c2_gmsl & valid_ucb_gmsl)) + mean_ucb ;
% gmsl3 = gmsl3(~isnan(gmsl3)); % 
% gmsl = gmsl-mean(gmsl(1:length(gmsl3)))+mean(gmsl3); % ????

% repeat for PODAAC data (which is GSFC)
yy = load([path2GMSL,'GMSL_GODDARD_2014_SEAS.txt']); % GMSL PO.DAAC
% see ftp://podaac.jpl.nasa.gov/allData/merged_alt/L2/TP_J1_OSTM/global_mean_sea_level/
gmsl3_goddard = yy(:,6);
% gmsl3_goddard11 = yy(:,11);

ty_goddard = yy(:,3);
ty_goddard = datenum(fix(ty_goddard),1,1) + rem(ty_goddard,fix(ty_goddard)).*365.24;
gmsl3_goddard = gmsl3_goddard./10; % convert from mm to cm
gmsl3_goddard = interp1(ty_goddard,gmsl3_goddard,dates);

valid_goddard_gmsl =  vec(~isnan(gmsl3_goddard)) ;

% gmsl3_goddard = gmsl3_goddard(~isnan(gmsl3_goddard));

% % % % % % % gmsl3_goddard11 = gmsl3_goddard11./10; % convert from mm to cm
% % % % % % % gmsl3_goddard11 = interp1(ty_goddard,gmsl3_goddard11,dates);
% % % % % % % gmsl3_goddard11 = gmsl3_goddard11(~isnan(gmsl3_goddard11));
% % % % % % % valid_goddard11_gmsl =  vec(~isnan(gmsl3_goddard11)) ;

% change mean of GMSL from NOAA to have same mean as UCB obviously this
% may change if we continue to use NOAA rather than UCB
% change mean of GMSL from NOAA to have same mean as UCB
gmsl3_goddard = gmsl3_goddard - mean(gmsl3_goddard(valid_goddard_gmsl & valid_ucb_gmsl)) + mean_ucb ;
% % % % % gmsl3_goddard11 = gmsl3_goddard11 - mean(gmsl3_goddard11(valid_goddard11_gmsl & valid_ucb_gmsl)) + mean_ucb ;


%   save /noc/users/cryo/QCV_Cryo2/code/validation/code_valid/testdata.mat dates gmsl* ty*  % y v



dl = datestr(dates);
dl = cellstr(dl);
if length(dl)>5
    dt = round(length(dl)/6);
else
    dt = 1;
end
dl2 = dl(2:dt:end);
idl = 2:dt:length(dl);
kdl = cellfun(@(x) strcmp(x,dl2{end}),dl);
kdl = find(kdl);
if (length(dl)-kdl) < dt/2 && length(dl) > 5;
    dl2(end) = dl(end);
    idl(end) = length(dl);
elseif length(dl) > 5
    dl2(end+1) = dl(end);
    idl(end+1) = length(dl);
else
    idl = 1:length(dl);
    dl2 = dl;
end
disp(gmsl(end))
    
rgb1 = [215/255 75/255 75/255];
rgb2 = [120/255 121/255 116/255];
rgb3 = [50/255 121/255 255/255];

hold on
h1 = plot(gmsl3,'-o','LineWidth',2,'MarkerFaceColor',rgb1,'color',rgb1);
h2 = plot(gmsl,'-o','lineWidth',2,'color',rgb2,'MarkerFaceColor',rgb2);
h3 = plot(gmsl3_goddard,'-o','lineWidth',2,'color',rgb3,'MarkerFaceColor',rgb3);
% % % % h4 = plot(gmsl3_goddard11,':+','lineWidth',2,'color',rgb3,'MarkerFaceColor',rgb3);

hl = legend([h1,h2 h3 ],'Colorado','CryoSat-2','GSFC','location','northWest');
set(gca,'FontSize',16)
ylabel('GMSL (cm)')
ylim([min(min(gmsl),min(gmsl3))-0.5 max(max(gmsl),max(gmsl3))+0.5])
set(gca,'XTick',idl,'XTickLabel',dl2)
set(hl,'box','off')
box on
xlim([0 length(gmsl)+1])
rotateXLabels(gca,45)
set(gcf,'paperPositionMode','auto')
print([dir_out 'Fig_GMSL_C2_CU'],'-r300','-depsc')
close all

    
    

% % interpolation
% z = NaN(nmonths,length(lat(:)));
% for i=1:nmonths
%     %F = scatteredInterpolant(x2{i}',y2{i}',v2{i}','natural');
%     F = TriScatteredInterp(x2{i}',y2{i}',v2{i}','natural');
%     z(i,:) = F(lon(:),lat(:));
% end
% 
% % restrict coordinates to (-180,180)
% k = find(lon(:,1)>=-180 & lon(:,1)<=180);
% z = reshape(z,size(z,1),size(lon,1),size(lon,2));
% z = z(:,k,:);
% lon = lon(k,:);
% lat = lat(k,:);
% z = z(:,:);
% lon = lon(:)';
% lat = lat(:)';
% k65 = abs(lat) <= 65;
% lon = lon(k65); %#ok
% lat = lat(k65);
% z = z(:,k65);
% gmsl2 = zeros(nmonths,1);
% for i=1:nmonths
%     knan = z(i,:) == z(i,:);
%     W = cosd(lat(knan));
%     gmsl2(i) = sum(bsxfun(@times,z(i,knan),W))/sum(W);
% end

        


