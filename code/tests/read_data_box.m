do not think this is used in routine processing
clear all

dayDate = '20160302';
path1 = '/Volumes/scratch/general/cryosat/daily_stats/';
path2 = '/Volumes/noc/users/fmc1q07/QCV_Cryo2/code/tests/';

z = load([path1,'Stats_SIR_FDM_L2_',dayDate,'.mat']);

x = z.lon;
y = z.lat;
t = z.tutc;
sla = z.ssha;
swh = z.swh;
fa1 = z.validFlag_ssha;
fw1 = z.validFlag_swh;
fa2 = z.validNOC_ssha;
fw2 = z.validNOC_swh;

k1 = x >= -90 & x <= 30 & y >= 0 & y <= 60;
k2 = x >= 120 & x <= 180 & y >= 0 & y <= 60;

ka1 = k1 & fa1 & fa2;
ka2 = k2 & fa1 & fa2;
kw1 = k1 & fw1 & fw2;
kw2 = k2 & fw1 & fw2;

xa1 = x(ka1);
xa2 = x(ka2);
xw1 = x(kw1);
xw2 = x(kw2);

ya1 = y(ka1);
ya2 = y(ka2);
yw1 = y(kw1);
yw2 = y(kw2);

ta1 = (t(ka1) - datenum(2000,1,1,0,0,0)).*24*3600; % seconds since 00:00:00 UTC on 1 Jan 2000 
ta2 = (t(ka2) - datenum(2000,1,1,0,0,0)).*24*3600;
tw1 = (t(kw1) - datenum(2000,1,1,0,0,0)).*24*3600;
tw2 = (t(kw2) - datenum(2000,1,1,0,0,0)).*24*3600;

a1 = sla(ka1);
a2 = sla(ka2);
w1 = swh(kw1);
w2 = swh(kw2);

oa1 = [ta1' ya1' xa1' a1'];
oa2 = [ta2' ya2' xa2' a2'];
ow1 = [tw1' yw1' xw1' w1'];
ow2 = [tw2' yw2' xw2' w2'];

if exist(dayDate,'dir') ~= 7
    mkdir(dayDate)
end

save([path2,dayDate,'/','ssha_atlantic_',dayDate,'.txt'],'oa1','-ascii','-double')
save([path2,dayDate,'/','ssha_pacific_',dayDate,'.txt'],'oa2','-ascii','-double')
save([path2,dayDate,'/','swh_atlantic_',dayDate,'.txt'],'ow1','-ascii','-double')
save([path2,dayDate,'/','swh_pacific_',dayDate,'.txt'],'ow2','-ascii','-double')

if exist([dayDate,'.zip'],'file') == 2
    delete([dayDate,'.zip'])
end
zip(dayDate,dayDate)
