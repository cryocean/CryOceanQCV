

function plot_xml_version_mode_mask(data_in)
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


hp = patchm(yh,xh,-1e+9,[0.7 0.7 0.7]);
set(hp,'EdgeColor','none')
lakes = shaperead('worldlakes', 'UseGeoCoords', true);
geoshow(ha,lakes, 'FaceColor', [1 1 1],'EdgeColor',[1 1 1])
% ADD POLYGONS AND RECTANGLES


X = data_in.X ;
Y = data_in.Y ;

zonecentre = data_in.centre ;


ksar = find(strcmp(data_in.mode  , 'SAR') & ~strcmp(    data_in.shape ,  'CIRCLE'  )) ;
ksin = find(strcmp(data_in.mode  , 'SIN') & ~strcmp(    data_in.shape ,  'CIRCLE'  ) ) ;
klrm = find(strcmp(data_in.mode  , 'LRM')  & ~strcmp(    data_in.shape ,  'CIRCLE'  )) ;

ksar_circ = find(strcmp(data_in.mode  , 'SAR') & strcmp(    data_in.shape ,  'CIRCLE'  )) ;
ksin_circ = find(strcmp(data_in.mode  , 'SIN') & strcmp(    data_in.shape ,  'CIRCLE'  ) ) ;
klrm_circ = find(strcmp(data_in.mode  , 'LRM')  & strcmp(    data_in.shape ,  'CIRCLE'  )) ;



ksar(end) =[];
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


% add in poles

linem(Y{end},X{end},'Color','k', ...
    'LineStyle','-','LineWidth',2);



framem
mlabel('MLabelParallel','south','MLabelLocation',-180:60:180)
plabel('PLabelLocation',-90:30:90)
setm(ha,'FontSize',12)
geoshow('landareas.shp','FaceColor','none','EdgeColor','black')
tightmap









