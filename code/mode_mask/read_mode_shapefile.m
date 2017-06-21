% reads mode mask in from XML not shapefile !
FOR REAL VERSION THAT READS XML NOT SHAPEFILE SEE read_mode_xml

clear;clc;close all
cd /noc/users/cryo/QCV_Cryo2/code/mode_mask/

addpath ~/matlabfns/
% % % % % S = shaperead('mask3 _8.shp');
% % % % % S = shaperead('mask3_9.shp');

% fn ='Mask_3_8.xml';

fn =  'mask3_9.xml' ;
% % % % % % DOMnode = xmlread('mask3_9.xml') ;
s = xml2struct(fn );
regs = s.SatX.MapRegions.MapRegion  ;

numregs = length(regs); zoneid = cell(numregs,1);
zonelabel = zoneid; zonetrans = zoneid; zoneshape = zoneid;
zonepoints = zoneid; zonecentre = NaN(numregs,2);
zoneradius = NaN(numregs,1);
zonepriority = zoneradius;

Y = zoneid; X = zoneid;

for n=1:numregs ;
    zoneid{n} = regs{1,n}.Parameters.ZoneId.Text;
    if isfield(regs{1,n}.Parameters, 'Description')
        zonelabel{n} = regs{1,n}.Parameters.Description.Text;
    else
        zonelabel{n} = zoneid{n} ;
    end
    zonepriority(n) = str2double( regs{1,n}.Parameters.Priority.Text);
    zonetrans{n} = regs{1,n}.Parameters.Transition.Text;
    zoneshape{n} = regs{1,n}.Shape.Text;
    zonecentre(n,1) = str2double(regs{1,n}.Center.Attributes.Lat);
    zonecentre(n,2) = str2double(regs{1,n}.Center.Attributes.Lon);
    if strcmp( 'CIRCLE',zoneshape(n) );
        zoneradius(n) = str2double(regs{1,n}.CircleRadius.Text);
    else
        zonepoints{n} = str2num(regs{1,n}.MapPoints.Text);
        Y{n} = zonepoints{n}(:,1) ;
        X{n} = zonepoints{n}(:,2) ;
    end
end; clear n

% % % % % % % % % % % % % % % % % % % % % % KML version
% % % % % % % % % % % % % % % % % % % % % addpath ~/matlabfns/
% % % % % % % % % % % % % % % % % % % % % s_xml = kml2struct( '/noc/users/cryo/QCV_Cryo2/code/mode_mask/Mask_3_9.kml' );
% % % % % % % % % % % % % % % % % % % % % warning('KML only has polygons and no circles?!??!!?')

% find LRM etc
isLRM = strcmp(zonetrans,'TO_LRM__') ;
isSAR = strcmp(zonetrans,'TO_SAR__') ;
isSIN = strcmp(zonetrans,'TO_SIN__') ;
% find shape type
isPOLYGON = strcmp(zoneshape,'POLYGON') ;
isCIRCLE = strcmp(zoneshape,'CIRCLE') ;
isRECTANGLE = strcmp(zoneshape,'RECTANGLE') ;

% find where shape and mode match (rectangle considered special case of
% polygon)
klrm = find(isLRM & (isPOLYGON | isRECTANGLE)) ;
ksar = find(isSAR & (isPOLYGON | isRECTANGLE));
ksin = find(isSIN & (isPOLYGON | isRECTANGLE));

klrm_circ = find( isLRM & isCIRCLE) ;
ksar_circ = find(isSAR & isCIRCLE );
ksin_circ = find(isSIN & isCIRCLE );

% checks have everything balanced
if (sum(isSAR) + sum(isLRM) + sum(isSIN)) ~= numregs ,error('types do not sum'), end;
if (sum(isPOLYGON) + sum(isCIRCLE) + sum(isRECTANGLE)) ~= numregs ,error('shapes do not sum'), end;
if (sum(length(klrm)) + sum(length(ksar)) + sum(length(ksin)) + ...
        sum(length(klrm_circ)) + sum(length(ksar_circ)) + ...
        sum(length(ksin_circ))) ~= numregs ;
    error('shapes/types do not sum');
end;

% plot
path4 = '/noc/users/cryo/QCV_Cryo2/code/';
load([path4 'gshhs_c.mat'])
xh = gshhs_c.lon;
yh = gshhs_c.lat;
lev = gshhs_c.level;
kl = (lev == 1 | lev == 6);
xh = xh(kl);
yh = yh(kl);
xh = cell2mat(xh);
yh = cell2mat(yh);
xh = xh(:)';
yh = yh(:)';

figure
pPos = [16.8 11.04];
pos = [12.6 10.8];
% maximise_window_cb
set(gcf,'PaperUnits', 'centimeters','PaperPosition', [0.5 0.5 pPos]);
set(gca,'Units','centimeters','Position',[1.1 0 pos])
ha = axesm('MapProjection','gstereo');
hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])
% ADD POLYGONS AND RECTANGLES
hold on
for i=1:length(ksar)
    h1 = patchm(Y{ksar(i)},X{ksar(i)},'Facecolor',[1 197/255 108/255], ...
        'EdgeColor','none','LineStyle','none');
end; clear i

hold on
for i=1:length(ksin)
    h2 = patchm(Y{ksin(i)},X{ksin(i)},'Facecolor',[110/255 197/255 233/255], ...
        'EdgeColor','none','LineStyle','none');
end; clear i

hold on
% for i=1:length(ksar2)
%     h1 = patchm(Y{ksar2(i)},X{ksar2(i)},'Facecolor',[1 197/255 108/255], ...
%         'EdgeColor','none','LineStyle','none'); %#ok
% end
%
% hold on
for i=1:length(klrm)
    h3 = patchm(Y{klrm(i)},X{klrm(i)},'Facecolor',[1 89/255 89/255], ...
        'EdgeColor','none','LineStyle','none');
end; clear i
% ADD CIRCLES
hold on
for i=1:length(ksar_circ)
    plotm(zonecentre(ksar_circ(i),:) , 'Marker','o',...
        'Color',[1 197/255 108/255],'MarkerSize',10,...
        'MarkerFaceColor',[1 197/255 108/255])
    plotm(zonecentre(ksar_circ(i),:) , 'Marker','*',...
        'Color','k','MarkerSize',12)
end ; clear i

hold on
for i=1:length(ksin_circ)
    plotm(zonecentre(ksin_circ(i),:) , 'Marker','o',...
        'Color',[110/255 197/255 233/255],'MarkerSize',10,...
        'MarkerFaceColor',[110/255 197/255 233/255])
    plotm(zonecentre(ksin_circ(i),:) , 'Marker','*',...
        'Color','k','MarkerSize',12)
end; clear i

hold on
for i=1:length(klrm_circ)
    plotm(zonecentre(klrm_circ(i),:) , 'Marker','o',...
        'Color', [1 89/255 89/255],'MarkerSize',10,...
        'MarkerFaceColor',[1 89/255 89/255])
    plotm(zonecentre(klrm_circ(i),:) , 'Marker','*',...
        'Color','k','MarkerSize',12)
end; clear i

framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
geoshow('landareas.shp','FaceColor','none','EdgeColor','black')
tightmap

%create baseline
clear gshhs_c h1 h2 h3
% recode zonetrans to mode
% SAR
zonetrans(strcmp(  zonetrans , 'TO_SAR__')) = {'SAR'} ;
zonetrans(strcmp(  zonetrans , 'TO_SIN__')) = {'SIN'} ;
zonetrans(strcmp(  zonetrans , 'TO_LRM__')) = {'LRM'} ;
% field1 = 'X'; field2 = 'Y'; field3 = 'mode'; field4 = 'priority';
% value1 = X ; value2 = Y; value3 = zonetrans; value4 = zonepriority;
% clear X Y zonetrans zonepriority
% struc_baseline = struct('X' ,X, 'Y' ,Y, 'mode' ,zonetrans, 'priority' , ...
%     zonepriority) ;

struc_baseline =  struct ;
struc_baseline.X = X;
struc_baseline.Y = Y;
struc_baseline.mode = zonetrans';
struc_baseline.priority = zonepriority';

% for each of the 24 seasonal variations of sea ice
% seasMask = struct ;
% seasMaskModes = struct('X' ,{}, 'Y' ,{}, 'mode' ,{}, 'priority' , ...
%     []) ;

% identify which are the poles
isArc = find(strcmp(zoneid,'CYSSAR02')) ;
isAnt = find(strcmp(zoneid,'CYSSAR01')) ;

struc_baseline.X([isArc isAnt]) = [];
struc_baseline.Y([isArc isAnt]) = [];
struc_baseline.mode([isArc isAnt]) = [];
struc_baseline.priority([isArc isAnt]) = [];

putPolesHere = length(struc_baseline.X) + 1;
% struc_baseline.X(putPolesHere) = {[]};
% struc_baseline.Y(putPolesHere) = {[]};
% struc_baseline.mode(putPolesHere) = {[]};
% struc_baseline.priority(putPolesHere) = 0;



for nt = 1:24 ;
    eval(['seasMaskModes.seas_'   sprintf('%02d',nt)  ' = struc_baseline ;']);
    % need to read in data for Arctic and Antarctic from v3_8 of mask (ASSUMES does not change wrt year)
    % seasMaskModes.seas01(put
    data_out = rtn_month_artant(nt) ;
    
    eval(['seasMaskModes.seas_' sprintf('%02d',nt) '.X{putPolesHere} = data_out(1,:)'' ;'])
    eval(['seasMaskModes.seas_' sprintf('%02d',nt) '.Y{putPolesHere} = data_out(2,:)'' ;'])
    eval(['seasMaskModes.seas_' sprintf('%02d',nt) '.mode{putPolesHere} = ''SAR'' ;'])
    eval(['seasMaskModes.seas_' sprintf('%02d',nt) '.priority(putPolesHere) = 99 ;'])
    
    clear data_out
    
%     figure
%     title(['nt is ' num2str(nt)])
%     plot(seasMaskModes.seas_01.X{putPolesHere}, seasMaskModes.seas_01.Y{putPolesHere},'k.')
    
end; clear nt










