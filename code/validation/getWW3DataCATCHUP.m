function getWW3DataCATCHUP
% get data from SOEST:
% (http://oos.soest.hawaii.edu/erddap/griddap/NWW3_Global_Best.html)
% UPDATED JAN 2017
path_ww3 = '/noc/mpoc/cryo/cryosat/validation_data_2/WWIII/'; % path where data are saved

% path_ww3 = '/scratch/general/cryosat/validation_data_2/WWIII/'; % path where data are saved
% list files with Matlab data in above folder
% fava = dir([path_ww3,'swh*.mat']);
% fava = struct2cell(fava);
% fava = fava(1,:);
% fava = cellfun(@(x) x(10:15),fava,'unif',0); 
% tNum = datenum(fava,'yyyymm') + 15; % convert from 1st of month to mid-point (approx) - 15th
% tfirst = max(tNum) + 30; % find month of first data still required


tfirst = datenum('20130701','yyyymmdd') ;
tlast = datenum('20130930','yyyymmdd') ;


% tlast = datenum(date,'dd-mmm-yyyy')-14; % latency of about 12 days in the data (from today's date)
% tlast = datestr(tlast,'mmm-yyyy');
% tlast = datenum(tlast,'mmm-yyyy')-1;

if tlast > tfirst
    months = tfirst:tlast; % need to download these months
    months = cellstr(datestr(months,'mmm-yyyy'));
    months = unique(months);
    for i=1:length(months)
        [Y,M] = datevec(datenum(months{i},'mmm-yyyy'));
        df = num2str(eomday(Y,M));
        Y = num2str(Y);
        M = num2str(M,'%02i');
        url = ['http://oos.soest.hawaii.edu/erddap/griddap/NWW3_Global_Best.mat?Thgt',...
            '[(' Y '-' M '-01T00:00:00Z):1:(' Y '-' M '-' df 'T23:59:59Z)]',...
            '[(0.0):1:(0.0)][(-77.5):1:(77.5)][(0.0):1:(359.5)]'];
        urlwrite(url,[path_ww3 'swh_nww3_' Y M '.mat']);
    end
end