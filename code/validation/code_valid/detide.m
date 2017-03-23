function [res,tide] = detide(y,t)
addpath('/noc/users/cryo/QCV_Cryo2/code/validation/code_valid/t_tide_v1')

y = y(:);
t = t(:);
[year,~,~] = datevec(t);
knan = (y == y);
yearU = unique(year(knan));

ny = numel(yearU);
tide = NaN(length(y),1);
res = NaN(length(y),1);
i = 1;
while i<=ny
    ku = find(year == yearU(i));
    ytmp = y(ku);
    if sum(ytmp == ytmp) < 4000
        if i ~= ny && i ~= 1
            ktmp = [ku;find(year == yearU(i+1))]; 
            ytmp = y(ktmp);
            if sum(ytmp == ytmp) < 4000
                ktmp = [find(year == yearU(i-1));ku];
                ytmp = y(ktmp);
            end
            ku = ktmp;
        elseif i == 1
            ktmp = [ku;find(year == yearU(i+1))]; 
            ytmp = y(ktmp);
            ku = ktmp;
        elseif i == ny
            ktmp = [find(year == yearU(i-1));ku]; 
            ytmp = y(ktmp);
            ku = ktmp;
        end
        if sum(ytmp == ytmp) < 4000
            i = i+1;
            continue
        end
    end
    k0 = find(ytmp == ytmp,1,'first');
    kf = find(ytmp == ytmp,1,'last');
    ytmp = ytmp(k0:kf);
    ku = ku(k0:kf);
    i = i+1;
    
    t0 = t(ku(1));
    nl = length(ytmp);
    mu = nanmean(ytmp);
    [TIDESTRUC,~] = t_tide(ytmp,'synthesis',3,'start time',t0, ...
        'output','none');
    T = t0+(0:1/24:(nl/24-1/24));
    xout = t_predic(T,TIDESTRUC,'synthesis',3);
    tide(ku) = xout+mu;
    res(ku) = ytmp-xout';
end
