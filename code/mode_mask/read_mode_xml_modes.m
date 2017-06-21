% reads mode mask in from XML not subscript of maskFromKML_generic_func.m
% that reads in and saves the required names etc
function modes = read_mode_xml_modes(version)





fn = ['mask' version '.xml'];


s = xml2struct(fn );


regs = s.SatX.MapRegions.MapRegion  ;

numregs = length(regs); zoneid = cell(numregs,1);
zonelabel = zoneid; zonetrans = zoneid; zoneshape = zoneid;
% zonepoints = zoneid; zonecentre = NaN(numregs,2);
% zoneradius = NaN(numregs,1);
zonepriority = NaN(numregs,1);

% Y = zoneid; X = zoneid;

for n=1:numregs ;
    zoneid{n} = regs{1,n}.Parameters.ZoneId.Text;
    if isfield(regs{1,n}.Parameters, 'Description')
        zonelabel{n} = regs{1,n}.Parameters.Description.Text;
    else
        zonelabel{n} = zoneid{n} ;
    end
    zonepriority(n) = str2double( regs{1,n}.Parameters.Priority.Text);
    zonetrans{n} = regs{1,n}.Parameters.Transition.Text;
    %     zoneshape{n} = regs{1,n}.Shape.Text;
    %     zonecentre(n,1) = str2double(regs{1,n}.Center.Attributes.Lat);
    %     zonecentre(n,2) = str2double(regs{1,n}.Center.Attributes.Lon);
    %     if strcmp( 'CIRCLE',zoneshape(n) );
    %         zoneradius(n) = str2double(regs{1,n}.CircleRadius.Text);
    %     else
    %         zonepoints{n} = str2num(regs{1,n}.MapPoints.Text);
    %         Y{n} = zonepoints{n}(:,1) ;
    %         X{n} = zonepoints{n}(:,2) ;
    %     end
end; clear n

% % % % % % % % % % % % % % % % % % % % % % KML version
% % % % % % % % % % % % % % % % % % % % % addpath ~/matlabfns/
% % % % % % % % % % % % % % % % % % % % % s_xml = kml2struct( '/noc/users/cryo/QCV_Cryo2/code/mode_mask/Mask_3_9.kml' );
% % % % % % % % % % % % % % % % % % % % % warning('KML only has polygons and no circles?!??!!?')

% find LRM etc
isLRM = strcmp(zonetrans,'TO_LRM__') ;
isSAR = strcmp(zonetrans,'TO_SAR__') ;
isSIN = strcmp(zonetrans,'TO_SIN__') ;
% % % % % % % % % % % % % % % % % % % % find shape type
% % % % % % % % % % % % % % % % % % % isPOLYGON = strcmp(zoneshape,'POLYGON') ;
% % % % % % % % % % % % % % % % % % % isCIRCLE = strcmp(zoneshape,'CIRCLE') ;
% % % % % % % % % % % % % % % % % % % isRECTANGLE = strcmp(zoneshape,'RECTANGLE') ;

% % % % % % % % % % % % % % % % % % % find where shape and mode match (rectangle considered special case of
% % % % % % % % % % % % % % % % % % % polygon)
% % % % % % % % % % % % % % % % % % klrm = find(isLRM & (isPOLYGON | isRECTANGLE)) ;
% % % % % % % % % % % % % % % % % % ksar = find(isSAR & (isPOLYGON | isRECTANGLE));
% % % % % % % % % % % % % % % % % % ksin = find(isSIN & (isPOLYGON | isRECTANGLE));
% % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % klrm_circ = find( isLRM & isCIRCLE) ;
% % % % % % % % % % % % % % % % % % ksar_circ = find(isSAR & isCIRCLE );
% % % % % % % % % % % % % % % % % % ksin_circ = find(isSIN & isCIRCLE );

% checks have everything balanced
if (sum(isSAR) + sum(isLRM) + sum(isSIN)) ~= numregs ,error('types do not sum'), end;
% % % % % % % % % % % % % % % % % % % % % % % if (sum(isPOLYGON) + sum(isCIRCLE) + sum(isRECTANGLE)) ~= numregs ,error('shapes do not sum'), end;
% % % % % % % % % % % % % % % % % % % % % % % if (sum(length(klrm)) + sum(length(ksar)) + sum(length(ksin)) + ...
% % % % % % % % % % % % % % % % % % % % % % %         sum(length(klrm_circ)) + sum(length(ksar_circ)) + ...
% % % % % % % % % % % % % % % % % % % % % % %         sum(length(ksin_circ))) ~= numregs ;
% % % % % % % % % % % % % % % % % % % % % % %     error('shapes/types do not sum');
% % % % % % % % % % % % % % % % % % % % % % % end;

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

modes2 =  cell(length(zonetrans),3) ;
% % % % % % % % % % % % struc_baseline.X = X;
% % % % % % % % % % % % struc_baseline.Y = Y;
% modes2 = zonetrans';
% modes2.priority = zonepriority';
% modes2.name = zoneid ;
% % % % % % % % % % struc_baseline.shape = zoneshape' ;
% % % % % % % % % % struc_baseline.centre(:,1) = zonecentre(:,1) ;
% % % % % % % % % % struc_baseline.centre(:,2) = zonecentre(:,2) ;
% % % % % % % % % % struc_baseline.radius = zoneradius';

% for each of the 24 seasonal variations of sea ice
% seasMask = struct ;
% seasMaskModes = struct('X' ,{}, 'Y' ,{}, 'mode' ,{}, 'priority' , ...
%     []) ;

for n = 1:length(zonetrans);
    
    
    modes2{n,1} = zoneid(n) ;
        modes2{n,2} = zonetrans(n) ;
        modes2{n,3} = zonepriority(n) ;

    
end; clear n

modes = cell2struct(modes2,{'name' 'mode' 'priority'},2) ;


