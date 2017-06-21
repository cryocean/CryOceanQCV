% function to return Arctic/Antarctic mode masks for given time period (nt)


function data_out = rtn_month_artant(nt)

% load mask v3_8

load /noc/users/cryo/QCV_Cryo2/code/mode_mask/seasMask_polar

data_out = polar{nt,1};



