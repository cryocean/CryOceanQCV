% script to produce new version of geographical mode mask


clear;close all;clc

cd /noc/users/cryo/QCV_Cryo2/code/mode_mask/
versionmask = '3_4' ; % set version mask BE AWARE that the naming 
% convention varies and some files might need titlte changing or 
% copyying eg Mask versus mask _3_n versus 3_n, etc


% % % % % % % % % 
read_mode_generic_fn(versionmask)
% % % % % % % % % 
% % % % % % % % % 
% % % % % % % % % 
% % % % % % % % % work on from here
% % % % % % % % % 
% % % % % % % % % 

load(['seasMaskAllModes_',versionmask,'.mat']);

for seas = 1:24;
    % read that mask then divide data into matrices
    mask = seasMaskModes.(['seas_' sprintf('%1.2d',seas)]);
    Mode = mask.mode;
    X = mask.X;
    Y = mask.Y;
    p = mask.priority;
    clear mask
    ksar = cellfun(@strcmp,Mode,repmat({'SAR'},1,length(X)));
    ksin = cellfun(@strcmp,Mode,repmat({'SIN'},1,length(X)));
    klrm = cellfun(@strcmp,Mode,repmat({'LRM'},1,length(X)));
    kp = p > 6 & ksar;
    
    xin = X(ksin);
    yin = Y(ksin);
    xsar = X(ksar & ~kp);
    ysar = Y(ksar & ~kp);
    xlrm = X(klrm);
    ylrm = Y(klrm);
    xsar2 = X(kp);
    ysar2 = Y(kp);
    
    % compute unions betwee SARin polygons
    K = cell(length(xin),1);
    for i=1:length(xin)
        for j=1:length(xin)
            [xi,~] = polybool('&',xin{i},yin{i},xin{j},yin{j});
            if ~isempty(xi)
                K{i} = [K{i} j];
            end
        end
    end
    
    xi = [];
    yi = [];
    for i=1:length(xin)
        for j=1:length(K{i})
            [xi,yi] = polybool('union',xi,yi,xin(K{i}(j)),yin(K{i}(j)));
        end
    end
    [xin,yin] = poly2cw(xi,yi);
    
    
    % compute unions between SAR polygons
    K = cell(length(xsar),1);
    for i=1:length(xsar)
        for j=1:length(xsar)
            [xi,~] = polybool('&',xsar{i},ysar{i},xsar{j},ysar{j});
            if ~isempty(xi)
                K{i} = [K{i} j];
            end
        end
    end
    
    xi = [];
    yi = [];
    for i=1:length(xsar)
        for j=1:length(K{i})
            [xi,yi] = polybool('union',xi,yi,xsar(K{i}(j)),ysar(K{i}(j)));
        end
    end
    [xsar,ysar] = poly2cw(xi,yi);
    
    % remove LRM from SAR
    for i=1:length(xsar)
        for j=1:length(xlrm)
            [xsar{i},ysar{i}] = polybool('minus',xsar{i},ysar{i},xlrm{j},ylrm{j});
        end
    end
    
    % remove LRM from SARin
    for i=1:length(xin)
        for j=1:length(xlrm)
            [xin{i},yin{i}] = polybool('minus',xin{i},yin{i},xlrm{j},ylrm{j});
        end
    end
    
    % remove SAR2 from SARin
    for i=1:length(xin)
        for j=1:length(xsar2)
            [xin{i},yin{i}] = polybool('minus',xin{i},yin{i},xsar2{j},ysar2{j});
        end
    end
    [yin,xin] = polyjoin(yin,xin);
    
    xin(end+1) = NaN; %#ok we need to add NaN so the isLRM function performs as expected
    yin(end+1) = NaN; %#ok
    z = [xin';yin'];
    SARinMask.(['seas_' sprintf('%1.2d',seas)]) = z;
end
save(['Mask_SARin_',versionmask,'.mat'],'SARinMask','-mat')
