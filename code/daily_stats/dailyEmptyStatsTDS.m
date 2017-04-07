function  dailyEmptyStatsTDS(day_date,data_type,outPath)

NOT NEEDED
fn = dir([outPath 'Stats_' data_type '*.mat']);



fldn = load([outPath fn(1).name]);
fldn = fieldnames(fldn);
dataStat = cell(length(fldn),1);
dataStat = cell2struct(dataStat,fldn);%#ok
% date_report = datestr(datenum(day_date,'yyyymmdd'),'dd/mm/yyyy');
% dataStat.date = {date_report};
save([outPath 'Stats_' data_type '_' day_date '.mat'],'-struct','dataStat')





end

