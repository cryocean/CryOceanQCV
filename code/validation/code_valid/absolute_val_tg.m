% absolute validation against tide gauges
function absolute_val_tg(month)
% input : month: string in format 'mmm-yyyy' ('Mar-2015') representing the
% year up to which the validation is performed

maxDist = 60; % maximum distance between altimetry obs and TG
dir_out = '/noc/users/cryo/QCV_Cryo2/code/gen_monthly_report/figures/';
% UPDATE HERE JAN 2017
% path_C2Data = '/scratch/general/cryosat/validation_data/';
path_C2Data = '/noc/mpoc/cryo/cryosat/validation_data/';
data_type = 'SIR_GOP_L2';

daten = datenum(month,'mmm-yyyy');
[Y,M] = datevec(daten);
Df = datestr([Y,M,eomday(Y,M),0,0,0],'yyyymmdd');

D0 = '20140410';
id_abs = [830 335 824 211 47 755 25]; %tgs used in absolute validation
tgLabel = {'La Coru�a(1)','Spring Bay(2)','Marseille(3)','Ponta Delgada(4)', ...
           'Chichijima(5)','Virginia Key(6)','Funafuti(7)'};
hellip = [51.9356 -4.6405 48.6042 55.7218 48.2670 -30.9693 33.6571]; % ellipsoidal heights
      
d0 = datenum(D0,'yyyymmdd');
df = datenum(Df,'yyyymmdd');
d = datestr(d0:df,'yyyymmdd');  
ndays = size(d,1);
files = [repmat(['dataVal_' data_type '_'],ndays,1) d,...
         repmat('.mat',ndays,1)];
      

ntg = length(id_abs);
c2Bias = cell(ntg,1);
for j=1:ntg
    disp(j)
    % read tide gauge data
    [xtg,ytg,ztg,ttg] = read_tgData_ascii(id_abs(j));
    
    [ztg,~] = detide(ztg,ttg); % remove tide
    ztg = ztg./1000; % convert to m
    
    %convert from tide-free system (GPS) to mean-tide system (altimetry)
    ctide = -0.12*(3/2*sind(ytg)^2-0.5);
    
    %convert to ellipsoidal height
    ztg = ztg + hellip(j) + ctide;
    
    %MSS at the tide gauge
    mssTg = find_mss(xtg,ytg);
    
    
    % read and select altimetry data
    x = NaN(ndays*1000,1);
    y = NaN(ndays*1000,1);
    v = NaN(ndays*1000,1);
    t = NaN(ndays*1000,1);
    o = NaN(ndays*1000,1);
    kit = 1;
    for i=1:ndays
        dataVal = load([path_C2Data files(i,:)]);
        d = Distance(xtg,ytg,dataVal.x_ssha,dataVal.y_ssha)./1000;
        k = d < maxDist;
        x(kit:(kit+sum(k)-1)) = dataVal.x_ssha(k);
        y(kit:(kit+sum(k)-1)) = dataVal.y_ssha(k);
        v(kit:(kit+sum(k)-1)) = dataVal.ssh(k) - dataVal.totOcn_tide(k) -...
                                dataVal.poleTide(k); %tides out, atm in
        t(kit:(kit+sum(k)-1)) = dataVal.t_ssha(k);
        o(kit:(kit+sum(k)-1)) = dataVal.orbit_ssha(k);
        kit = kit + sum(k);
    end
    k = find(~isnan(x));
    x = x(k);
    y = y(k);
    v = v(k);
    t = t(k);
    o = o(k);
    
    % temporal interpolation for the tg
    vtg = zeros(length(k),1);
    tb = zeros(length(k),1);
    for i=1:length(k)
        dt = abs(t(i)-ttg);
        kt = find(dt == min(dt),1);
        if(dt(kt) < 1)
            k0 = max(kt-10,1);
            kf = min(kt+10,length(ztg));
            vtg(i) = interp1(ttg(k0:kf),ztg(k0:kf),t(i));
            tb(i) = t(i);
        end
    end
    vtg(vtg == 0) = NaN;
    tb(vtg == 0 )= NaN;
    
    %hg = zeros(length(k),1);
    hbias = zeros(length(k),1);
    vv = zeros(length(k),1);
    for i=1:length(k)
        hmss = find_mss(x(i),y(i));
        hbias(i) = v(i)-vtg(i)-(hmss-mssTg);
        vv(i) = v(i) - (hmss-mssTg);
    end
    
    % select one value per pass
    ou = unique(o);
    hbias2 = NaN(size(ou));
    tb2 = NaN(size(ou));
    for i=1:length(ou)
        ko = find(o == ou(i) & ~isnan(hbias));
        if ~isempty(ko)
            kmin = find(abs(hbias(ko)) == min(abs(hbias(ko))),1);
            hbias2(i) = hbias(ko(kmin));
            tb2(i) = tb(ko(kmin));
        end
    end
    
    c2Bias{j} = hbias2;
    
end

mu = cellfun(@nanmean,c2Bias).*100;
STD = cellfun(@nanstd,c2Bias).*100;
muT = mean(mu);

ymax = max(mu+STD);
ymin = min(mu-STD);


% plot bias
figure(1)
pPos = [15 15];
pos = [10 6];
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.5 2.6 pos])
errorbar(mu,STD,'o','LineWidth',1.5,'color','k')
hold on
plot(0:8,muT.*ones(9,1),'LineWidth',2,'color','k')
plot(0:8,zeros(9,1),'--','LineWidth',1,'color','k')
text(4.5,ymax+1.5,['mean bias = ',num2str(muT,'%3.1f'),' cm'],'FontSize',12)
box on
xlim([0 8])
ylim([ymin-2,ymax+3])
ylabel('CryoSat-2 bias (cm)')
set(gca,'FontSize',12,'XTick',1:7,'XTickLabel',tgLabel,'YTick',-50:5:50)
rotateXLabels(gca(),45)

set(gcf, 'PaperPositionMode', 'manual')
print('-depsc',[dir_out 'Fig_c2_bias_tg'],'-r300') % save figure
close all






