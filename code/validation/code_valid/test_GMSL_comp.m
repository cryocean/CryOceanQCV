% short script to test algorithm for GMSL_computation, using first data
% in Sept 2016

clear;close all; clc
addpath /noc/mpoc/smos/software/tools_cb
load testdata

path2GMSL = '/noc/mpoc/cryo/cryosat/validation_data_2/';

lat2 = -65:1:65; % set latitude range
 ai = cosd(lat2);
lat_bin = NaN(size(lat2)) ;
num_cells = length(y);
xbar = NaN(num_cells,1) ;
sd_out =xbar ;
wmed = xbar ;


for n = 1:num_cells;
    ai_ind = cosd(y{n}) ;
    
     [xbar(n) ,sd_out(n), wmed(n)] = weighted_mean_etc_fn(...
         v{n},ai_ind)  ;

    
end; clear n


xbar = 100 .* xbar ;

% 
% for j=1:length(lat2)
%     k = save_lat >= lat2(j)-0.5 & save_lat<lat2(j)+0.5;
%     
%     lat_bin(k) = j ;
%     
% end
% 
% appears to differ
% 



cb_chart_set_size(21)

plot(gmsl,'k-+')
hold on
plot(xbar,'r-d')
maximise_window_cb
plot(gmsl3,'-mo')
plot(gmsl3_goddard11,'gv-')

legend('current','indi','UCB','GSFC')

cb_chart_set_size(21)


