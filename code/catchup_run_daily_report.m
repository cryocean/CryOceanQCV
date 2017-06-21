% file to run and test monthly_report_IOP_FDM_v3_nodisplay_CATCHUP
clear ;close all;clc
data_type='SIR_GOP_L2';


cd /noc/users/cryo/QCV_Cryo2/code
% addpath read_esa
warning('filenames of output have changed')

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


in_dir = '/noc/mpoc/cryo/cryosat/catchup/data_out/' ;
% out_dir  = '/noc/mpoc/cryo/cryosat/catchup/data_out/' ;

% % % % start = '20101201';
% % % % stop = '20140430';

start = '20130908';
stop = '20130910';

start_num = datenum(start , 'yyyymmdd') ;
stop_num = datenum(stop , 'yyyymmdd') ;


[YYMMDD(:,1),YYMMDD(:,2),YYMMDD(:,3)] = datevec(start_num:stop_num)  ;
months_use = unique(YYMMDD(:,1:2),'rows'); clear YYMMDD


for  nn = 1:size(months_use,1)  ;
    
    monthly_report_IOP_FDM_v3_nodisplay_CATCHUP(data_type,...
        datestr(datenum(months_use(nn,1),months_use(nn,2),1),'mmm-yyyy'))
    
end ; clear date_use;