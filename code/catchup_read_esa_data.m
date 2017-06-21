% file to run and test daily_stats_v6_FOR_CATCHUP
clear ;close all;clc
data_type='SIR_GOP_L2';


cd /noc/users/cryo/QCV_Cryo2/code
addpath read_esa


% % % % % % 
% % % % % % path1 =  ['/noc/mpoc/cryo/cryosat/catchup/' data_type '/' day_date(1:4) '/' day_date(5:6) '/'] ;%: path to the matlab files containing data to process
% % % % % % path2 = ['/scratch/general/cryosat/' data_type '/'];%: path to latency data
% % % % % % path3 = '/noc/mpoc/cryo/cryosat/daily_stats/' ;%;: path to folder where statistics are saved
% % % % % % mss  =  'DTU10MSS';%: DTU10 Mean Sea sSurfacee
% % % % % % gshhs = 'gshhs_i';%: coast lines (seems to be intermediate product)
% % % % % % ground = 'A';%: groundtracks
% % % % % % 
% % % % % % 
% % % % % % daily_stats_v6_FOR_CATCHUP(day_date,data_type,path1,path2,path3,mss,gshhs,ground)
% % % % % % 


in_dir = '/noc/mpoc/cryo/cryosat/catchup/SIR_GOP_L2/' ;
out_dir  = '/noc/mpoc/cryo/cryosat/catchup/data_out/' ;

start = '20101101';
stop = '20140430';



read_cryo2Data_catchup(in_dir,out_dir,start,stop)

