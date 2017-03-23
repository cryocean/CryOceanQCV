clear;close all;clc
cd /noc/users/cryo/QCV_Cryo2/code/


load for_test.mat

daily_stats_v6(t_window(i,:),data_type,path2,in_dir,path1, ...
    DTU10MSS,gshhs_i,A);
