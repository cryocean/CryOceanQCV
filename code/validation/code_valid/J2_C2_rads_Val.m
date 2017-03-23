% Comparison of Jason2 and CryoSat2
function J2_C2_rads_Val(month)
% input : month: string in format 'mmm-yyyy' ('Mar-2015')
% UPDATED JAN 2017
% path_C2Data = '/scratch/general/cryosat/validation_data/';
% path3 = '/scratch/general/cryosat/validation_data_2/'; %gshhs
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
path3 = '/noc/mpoc/cryo/cryosat/validation_data_2/'; %gshhs


pathOut = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/'; % for figures
path2 = '/noc/users/cryo/QCV_Cryo2/code/mode_mask/';
pathj2 = '/noc/mpoc/scratch/datasets/alt/rads/j2/a/';
pathc2 = '/noc/mpoc/scratch/datasets/alt/rads/c2/a/';
date0 = datenum(month,'mmm-yyyy');
[Y,M] = datevec(date0);
datef = date0+eomday(Y,M)-1;

days = datestr(date0:datef,'yyyymmdd');
ndays = size(days,1); % number of days in month


% ------------ polar regions based on the mode mask polygons --------------
load([path2 'seasMask_polar.mat']); % arctic mask

% choose seasonal mode mask
day0 = days(15,:);
seasMonth = day0(5:6);
seasDay = day0(end-1:end);
seas = 2*(str2num(seasMonth)-1); %#ok
if(str2double(seasDay) < 15)
    seas = seas+1;
else
    seas = seas+2;
end
polar = polar{seas}; %#ok

load([path3 'gshhs_i.mat']) % coastline
% -------------------------------------------------------------------------

slaJ2 = cell(ndays,1);
swhJ2 = cell(ndays,1);
wspJ2 = cell(ndays,1);
slaC2 = cell(ndays,1);
swhC2 = cell(ndays,1);
wspC2 = cell(ndays,1);
mhJ2 = 0;
mhC2 = 0;
msJ2 = 0;
msC2 = 0;
mwJ2 = 0;
mwC2 = 0;
nhJ2 = 0;
nhC2 = 0;
nsJ2 = 0;
nsC2 = 0;
nwJ2 = 0;
nwC2 = 0;
for nd=1:ndays
    disp(nd)
    day = days(nd,:);
    dirn = dir([pathj2 'c*']);
    fln = cell(10000,1); % preallocate large cell array then discard empty
    ct = 0;
    for i =1:length(dirn)
        fn = dir([pathj2 dirn(i).name '/*.nc']);
        ncid = netcdf.open([pathj2 dirn(i).name '/' fn(1).name],'NC_NOWRITE');
        time0 = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'first_meas_time');
        netcdf.close(ncid);
        ncid = netcdf.open([pathj2 dirn(i).name '/' fn(end).name],'NC_NOWRITE');
        timef = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'last_meas_time');
        netcdf.close(ncid)
        time0 = datenum(strrep(time0(1:10),'-',''),'yyyymmdd');
        timef = datenum(strrep(timef(1:10),'-',''),'yyyymmdd');
        if timef < datenum(day,'yyyymmdd')
            continue
        elseif time0 > datenum(day,'yyyymmdd')
            break
        else
            for j=1:length(fn)
                ncid = netcdf.open([pathj2 dirn(i).name '/' fn(j).name],'NC_NOWRITE');
                time0 = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'first_meas_time');
                timef = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'last_meas_time');
                netcdf.close(ncid)
                
                time0 = datenum(strrep(time0(1:10),'-',''),'yyyymmdd');
                timef = datenum(strrep(timef(1:10),'-',''),'yyyymmdd');
                if time0 <= datenum(day,'yyyymmdd') && timef >= datenum(day,'yyyymmdd')
                    ct = ct+1;
                    fln{ct} = [pathj2 dirn(i).name '/' fn(j).name];
                elseif time0 > datenum(day,'yyyymmdd')
                    break
                end
            end
        end
    end
    kempty = cellfun(@isempty,fln);
    fln(kempty) = [];
    
    SSHA = cell(length(fln),1);
    SWH = cell(length(fln),1);
    WSP = cell(length(fln),1);
    
    for i=1:length(fln)
        ncid = netcdf.open(fln{i},'NC_NOWRITE');
        t = netcdf.inqVarID(ncid,'time');
        t = netcdf.getVar(ncid,t,'double');
        
        % coordinates are used only to select ocen points
        lat = netcdf.inqVarID(ncid,'lat');
        fillVal = netcdf.getAtt(ncid,lat,'_FillValue','double');
        sf = netcdf.getAtt(ncid,lat,'scale_factor');
        lat = netcdf.getVar(ncid,lat,'double');
        lat(lat == fillVal) = NaN;
        lat = lat.*sf;
        
        lon = netcdf.inqVarID(ncid,'lon');
        fillVal = netcdf.getAtt(ncid,lon,'_FillValue','double');
        sf = netcdf.getAtt(ncid,lon,'scale_factor');
        lon = netcdf.getVar(ncid,lon,'double');
        lon(lon == fillVal) = NaN;
        lon = lon.*sf;
        lon(lon > 180) = lon(lon > 180) -360;
        % --------------------------------- coordinates
        
        ssha = netcdf.inqVarID(ncid,'ssha');
        fillVal = netcdf.getAtt(ncid,ssha,'_FillValue','double');
        sf = netcdf.getAtt(ncid,ssha,'scale_factor');
        ssha = netcdf.getVar(ncid,ssha,'double');
        ssha(ssha == fillVal) = NaN;
        ssha = ssha.*sf;
        
        
        swh = netcdf.inqVarID(ncid,'swh_ku');
        fillVal = netcdf.getAtt(ncid,swh,'_FillValue','double');
        sf = netcdf.getAtt(ncid,swh,'scale_factor');
        swh = netcdf.getVar(ncid,swh,'double');
        swh(swh == fillVal) = NaN;
        swh = swh.*sf;
        
        
        wsp = netcdf.inqVarID(ncid,'wind_speed_alt');
        fillVal = netcdf.getAtt(ncid,wsp,'_FillValue','double');
        sf = netcdf.getAtt(ncid,wsp,'scale_factor');
        wsp = netcdf.getVar(ncid,wsp,'double');
        wsp(wsp == fillVal) = NaN;
        wsp = wsp.*sf;
        
        netcdf.close(ncid);
        
        % ------- select points over ocean and outside polar regions ------
        knoPolar = isLRM(lon(:),lat(:),polar)';
        kocean = ~island(lon(:),lat(:),gshhs_i,1)' & abs(lat(:)) < 66;
        % --------------------------------------------------- end selecting
        
        t = datenum(1985,1,1)+t./(24*3600);
        kday = fix(t) == datenum(day,'yyyymmdd');
        SSHA{i} = ssha(kday & knoPolar & kocean);
        SWH{i} = swh(kday & knoPolar & kocean);
        WSP{i} = wsp(kday & knoPolar & kocean);
    end
    SSHA = cell2mat(SSHA);
    SWH = cell2mat(SWH);
    WSP = cell2mat(WSP);
    
    SSHA(abs(SSHA).*100 > 51) = NaN;
    [ch,xh] = hist(SSHA.*100,(-50:2:50));
    
    SWH(SWH > 12.125 | SWH < 0) = NaN;
    [cs,xs] = hist(SWH,(0:0.25:12));
    
    WSP(WSP > 24.25) = NaN;
    [cw,xw] = hist(WSP,(0:0.5:24));
    
    
    %     slaJ2{nd} = [ch;xh]; edited to be row vector 7/9/16
    slaJ2{nd} = [ch;xh(:)']; %
    %     swhJ2{nd} = [cs;xs]; edited to be row vector  7/2/17
    swhJ2{nd} = [cs;xs(:)'];
    %     wspJ2{nd} = [cw;xw]; edited to be row vector  7/2/17
    wspJ2{nd} = [cw;xw(:)'];
    
    mhJ2 = mhJ2 + sum(SSHA(~isnan(SSHA)));
    nhJ2 = nhJ2 + sum(~isnan(SSHA));
    msJ2 = msJ2 + sum(SWH(~isnan(SWH)));
    nsJ2 = nsJ2 + sum(~isnan(SWH));
    mwJ2 = mwJ2 + sum(WSP(~isnan(WSP)));
    nwJ2 = nwJ2 + sum(~isnan(WSP));
    
    
    % read CryoSat2 data
    datC2 = load([path_C2Data 'dataVal_SIR_GOP_L2_' day '.mat']);
    SSHA = datC2.ssha;
    SWH = datC2.swh;
    WSP = datC2.wsp;
    
    SSHA(abs(SSHA).*100 > 51) = NaN;
    [ch,xh] = hist(SSHA.*100,(-50:2:50));
    
    SWH(SWH > 12.125 | SWH < 0) = NaN;
    [cs,xs] = hist(SWH,(0:0.25:12));
    
    WSP(WSP > 24.25) = NaN;
    [cw,xw] = hist(WSP,(0:0.5:24));
    
    %     slaC2{nd} = [ch;xh]; %edited to be row vector 7/9/16
    slaC2{nd} = [ch;xh(:)'];
    swhC2{nd} = [cs;xs];
    wspC2{nd} = [cw;xw];
    
    mhC2 = mhC2 + sum(SSHA(~isnan(SSHA)));
    nhC2 = nhC2 + sum(~isnan(SSHA));
    msC2 = msC2 + sum(SWH(~isnan(SWH)));
    nsC2 = nsC2 + sum(~isnan(SWH));
    mwC2 = mwC2 + sum(WSP(~isnan(WSP)));
    nwC2 = nwC2 + sum(~isnan(WSP));
    
end
mhJ2 = mhJ2/nhJ2.*100;
mhC2 = mhC2/nhC2.*100;
msJ2 = msJ2/nsJ2;
msC2 = msC2/nsC2;
mwJ2 = mwJ2/nwJ2;
mwC2 = mwC2/nwC2;



% ----------------- do the same for CryoSat2 from RADS --------------------
slaC2r = cell(ndays,1);
swhC2r = cell(ndays,1);
wspC2r = cell(ndays,1);
mhC2r = 0;
msC2r = 0;
mwC2r = 0;
nhC2r = 0;
nsC2r = 0;
nwC2r = 0;
for nd=1:ndays
    disp(nd)
    day = days(nd,:);
    dirn = dir([pathc2 'c*']);
    fln = cell(10000,1); % preallocate large cell array then discard empty
    ct = 0;
    for i =1:length(dirn)
        fn = dir([pathc2 dirn(i).name '/*.nc']);
        ncid = netcdf.open([pathc2 dirn(i).name '/' fn(1).name],'NC_NOWRITE');
        time0 = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'first_meas_time');
        netcdf.close(ncid);
        ncid = netcdf.open([pathc2 dirn(i).name '/' fn(end).name],'NC_NOWRITE');
        timef = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'last_meas_time');
        netcdf.close(ncid)
        time0 = datenum(strrep(time0(1:10),'-',''),'yyyymmdd');
        timef = datenum(strrep(timef(1:10),'-',''),'yyyymmdd');
        if timef < datenum(day,'yyyymmdd')
            continue
        elseif time0 > datenum(day,'yyyymmdd')
            break
        else
            for j=1:length(fn)
                ncid = netcdf.open([pathc2 dirn(i).name '/' fn(j).name],'NC_NOWRITE');
                time0 = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'first_meas_time');
                timef = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'last_meas_time');
                netcdf.close(ncid)
                
                time0 = datenum(strrep(time0(1:10),'-',''),'yyyymmdd');
                timef = datenum(strrep(timef(1:10),'-',''),'yyyymmdd');
                if time0 <= datenum(day,'yyyymmdd') && timef >= datenum(day,'yyyymmdd')
                    ct = ct+1;
                    fln{ct} = [pathc2 dirn(i).name '/' fn(j).name];
                elseif time0 > datenum(day,'yyyymmdd')
                    break
                end
            end
        end
    end
    
    kempty = cellfun(@isempty,fln);
    fln(kempty) = [];
    
    SSHA = cell(length(fln),1);
    SWH = cell(length(fln),1);
    WSP = cell(length(fln),1);
    
    for i=1:length(fln)
        ncid = netcdf.open(fln{i},'NC_NOWRITE');
        t = netcdf.inqVarID(ncid,'time');
        t = netcdf.getVar(ncid,t,'double');
        
        % coordinates are used only to select ocen points
        lat = netcdf.inqVarID(ncid,'lat');
        fillVal = netcdf.getAtt(ncid,lat,'_FillValue','double');
        sf = netcdf.getAtt(ncid,lat,'scale_factor');
        lat = netcdf.getVar(ncid,lat,'double');
        lat(lat == fillVal) = NaN;
        lat = lat.*sf;
        
        lon = netcdf.inqVarID(ncid,'lon');
        fillVal = netcdf.getAtt(ncid,lon,'_FillValue','double');
        sf = netcdf.getAtt(ncid,lon,'scale_factor');
        lon = netcdf.getVar(ncid,lon,'double');
        lon(lon == fillVal) = NaN;
        lon = lon.*sf;
        lon(lon > 180) = lon(lon > 180) -360;
        % --------------------------------- coordinates
        
        ssha = netcdf.inqVarID(ncid,'ssha');
        fillVal = netcdf.getAtt(ncid,ssha,'_FillValue','double');
        sf = netcdf.getAtt(ncid,ssha,'scale_factor');
        ssha = netcdf.getVar(ncid,ssha,'double');
        ssha(ssha == fillVal) = NaN;
        ssha = ssha.*sf;
        
        swh = netcdf.inqVarID(ncid,'swh_ku');
        fillVal = netcdf.getAtt(ncid,swh,'_FillValue','double');
        sf = netcdf.getAtt(ncid,swh,'scale_factor');
        swh = netcdf.getVar(ncid,swh,'double');
        swh(swh == fillVal) = NaN;
        swh = swh.*sf;
        
        wsp = netcdf.inqVarID(ncid,'wind_speed_alt');
        fillVal = netcdf.getAtt(ncid,wsp,'_FillValue','double');
        sf = netcdf.getAtt(ncid,wsp,'scale_factor');
        wsp = netcdf.getVar(ncid,wsp,'double');
        wsp(wsp == fillVal) = NaN;
        wsp = wsp.*sf;
        
        netcdf.close(ncid);
        
        % ------- select points over ocean and outside polar regions ------
        knoPolar = isLRM(lon(:),lat(:),polar)';
        kocean = ~island(lon(:),lat(:),gshhs_i,1)' & abs(lat(:)) < 66;
        % --------------------------------------------------- end selecting
        
        t = datenum(1985,1,1)+t./(24*3600);
        kday = fix(t) == datenum(day,'yyyymmdd');
        SSHA{i} = ssha(kday & knoPolar & kocean);
        SWH{i} = swh(kday & knoPolar & kocean);
        WSP{i} = wsp(kday & knoPolar & kocean);
        
    end
    
    SSHA = cell2mat(SSHA);
    SWH = cell2mat(SWH);
    WSP = cell2mat(WSP);
    
    SSHA(abs(SSHA).*100 > 51) = NaN;
    [ch,xh] = hist(SSHA.*100,(-50:2:50));
    
    SWH(SWH > 12.125 | SWH < 0) = NaN;
    [cs,xs] = hist(SWH,(0:0.25:12));
    
    WSP(WSP > 24.25) = NaN;
    [cw,xw] = hist(WSP,(0:0.5:24));
    
    slaC2r{nd} = [ch(:)';xh(:)'];
    swhC2r{nd} = [cs(:)';xs(:)'];
    wspC2r{nd} = [cw(:)';xw(:)'];
    
    mhC2r = mhC2r + sum(SSHA(~isnan(SSHA)));
    nhC2r = nhC2r + sum(~isnan(SSHA));
    msC2r = msC2r + sum(SWH(~isnan(SWH)));
    nsC2r = nsC2r + sum(~isnan(SWH));
    mwC2r = mwC2r + sum(WSP(~isnan(WSP)));
    nwC2r = nwC2r + sum(~isnan(WSP));
    
end

mhC2r = mhC2r/nhC2r.*100;
msC2r = msC2r/nsC2r;
mwC2r = mwC2r/nwC2r;

% --------------------------------------------- end processing C2 RADS data




% ---------------------------- plot histograms ----------------------------
slaC2 = cell2mat(slaC2);
slaC2r = cell2mat(slaC2r);
slaJ2 = cell2mat(slaJ2);
swhC2 = cell2mat(swhC2);
swhC2r = cell2mat(swhC2r);
swhJ2 = cell2mat(swhJ2);
wspC2 = cell2mat(wspC2);
wspC2r = cell2mat(wspC2r);
wspJ2 = cell2mat(wspJ2);

[i1,~] = size(slaC2);
k = 1:2:i1;
xh = slaC2(2,:);
slaC2 = sum(slaC2(k,:),1);
slaC2r = sum(slaC2r(k,:),1);
slaJ2 = sum(slaJ2(k,:),1);

[i1,~] = size(swhC2);
k = 1:2:i1;
xs = swhC2(2,:);
swhC2 = sum(swhC2(k,:),1);
swhC2r = sum(swhC2r(k,:),1);
swhJ2 = sum(swhJ2(k,:),1);

[i1,~] = size(wspC2);
k = 1:2:i1;
xw = wspC2(2,:);
wspC2 = sum(wspC2(k,:),1);
wspC2r = sum(wspC2r(k,:),1);
wspJ2 = sum(wspJ2(k,:),1);

rgb2 = [51/255 153/255 1];
rgb1 = [215/255 75/255 75/255];
rgb3 = [156/255 207/255 49/255];

% -----------------SSHA
pPos = [15 13];
pos = [12 10.5];
figure(1)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[2.2 1.7 pos])
hold on
%hb1 = bar(x1,h1,'FaceColor',rgb2,'EdgeColor','none');
hb2 = bar(xh,slaC2/trapz(xh,slaC2),0.5,'FaceColor',rgb2,'EdgeColor','none');
hb1 = plot(xh,slaJ2/trapz(xh,slaJ2),'LineWidth',2.5,'color',rgb1);
hb3 = plot(xh,slaC2r/trapz(xh,slaC2r),'LineWidth',2.5,'color',rgb3);
box on
hl = legend([hb1 hb2 hb3],{['Jason-2 (\mu = ' num2str(mhJ2,'%4.2f') ' cm)'], ...
    ['GOP CryoSat-2 (\mu = ' num2str(mhC2,'%4.2f') ' cm)'],...
    ['RADS CryoSat-2 (\mu = ' num2str(mhC2r,'%4.2f') ' cm)']},'FontSize',16);
set(hl,'box','off','location','northwest')
set(gca,'FontSize',16,'XTick',-50:10:50)
xlim([-50 50])
xlabel('SSHA (cm)')
ylabel('Pdf')
yL = get(gca,'Ylim');
ylim([yL(1) 1.3*yL(2)])

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[pathOut 'Fig_hist_ssha_j2'],'-painters','-r300') % save figure
close all



% -----------------SWH
figure(1)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[2.2 1.7 pos])
hold on
%hb1 = bar(x1,h1,'FaceColor',rgb2,'EdgeColor','none');
hb2 = bar(xs,swhC2/trapz(xs,swhC2),0.5,'FaceColor',rgb2,'EdgeColor','none');
hb1 = plot(xs,swhJ2/trapz(xs,swhJ2),'LineWidth',2.5,'color',rgb1);
hb3 = plot(xs,swhC2r/trapz(xs,swhC2r),'LineWidth',2.5,'color',rgb3);
box on
hl = legend([hb1 hb2 hb3],{['Jason-2 (\mu = ' num2str(msJ2,'%4.2f') ' m)'], ...
    ['GOP CryoSat-2 (\mu = ' num2str(msC2,'%4.2f') ' m)'],...
    ['RADS CryoSat-2 (\mu = ' num2str(msC2r,'%4.2f') ' m)']},'FontSize',16);
set(hl,'box','off')
set(gca,'FontSize',16)
xlim([0 12])
xlabel('SWH (m)')
ylabel('Pdf')

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[pathOut 'Fig_hist_swh_j2'],'-painters','-r300') % save figure
close all


% -----------------WSP
figure(1)
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[2.2 1.7 pos])
hold on
%hb1 = bar(x1,h1,'FaceColor',rgb2,'EdgeColor','none');
hb2 = bar(xw,wspC2/trapz(xw,wspC2),0.5,'FaceColor',rgb2,'EdgeColor','none');
hb1 = plot(xw,wspJ2/trapz(xw,wspJ2),'LineWidth',2.5,'color',rgb1);
hb3 = plot(xw,wspC2r/trapz(xw,wspC2r),'LineWidth',2.5,'color',rgb3);

box on
hl = legend([hb1 hb2 hb3],{['Jason-2 (\mu = ' num2str(mwJ2,'%4.2f') ' m/s)'], ...
    ['GOP CryoSat-2 (\mu = ' num2str(mwC2,'%4.2f') ' m/s)'],...
    ['RADS CryoSat-2 (\mu = ' num2str(mwC2r,'%4.2f') ' m/s)']},'FontSize',16);
set(hl,'box','off')
set(gca,'FontSize',16)
xlim([0 24])
xlabel('Wind speed (m/s)')
ylabel('Pdf')
yL = get(gca,'Ylim');
ylim([yL(1) 1.35*yL(2)])

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[pathOut 'Fig_hist_wsp_j2'],'-painters','-r300') % save figure
close all
% -------------------------------------------------------------------------



