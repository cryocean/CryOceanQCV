% reads mode mask in from XML AND KML shapefile !
function read_mode_generic_fn(versionmask)

% cd /noc/users/cryo/QCV_Cryo2/code/mode_mask/
addpath /noc/users/cryo/matlabfns/ % xml and kml reading functions

% Read xml verison of Mode Mask into structure array s
fnx = dir(['*ask*' versionmask '.xml']);
fnx =  fnx(1).name;
warning('off','map:polygon:noExternalContours'); % prevents warnings that don't seem to affect results
s = xml2struct(fnx);
regs = s.SatX.MapRegions.MapRegion;
numregs = length(regs); % Save number of regions read

% Now read the kml version of the Mode Mask file for detailed polygons
fnk = dir(['*ask*' versionmask '.kml']);
fnk =  fnk(1).name;
s = kml_shapefile(fnk);

% Make sure we have the same number of regions
if length(s)~=numregs
  error(['xml file ' fnx ' has ' num2str(numregs) ' regions,'...
         'kml file ' fnk ' has ' num2str(length(s)) ' regions: check files!']);
end

% preallocate cell and double arrays arrays
% keep zoneid, priority and mode from xml file
Z = cell(numregs,3); % cell array used to convert to structure
zoneid = cell(numregs,1);
zonetrans = zoneid; zoneshape = zoneid; zoneprioritycell = zoneid;
% zonelabel = zoneid; zonepoints = zoneid;
% Y = zoneid; X = zoneid;

% zonecentre = NaN(numregs,2);
% zonepriority = NaN(numregs,1); zoneradius = zonepriority;

% Keep x and y from kml shapefile
x = zoneid; y = zoneid;

% transfer structure information to arrays
for n=1:numregs;
  % From xml file
    zoneid{n} = regs{1,n}.Parameters.ZoneId.Text;
%     if isfield(regs{1,n}.Parameters, 'Description')
%         zonelabel{n} = regs{1,n}.Parameters.Description.Text;
%     else
%         zonelabel{n} = zoneid{n} ;
%     end

    zoneprioritycell{n} = str2double(regs{1,n}.Parameters.Priority.Text);
    zonetrans{n} = regs{1,n}.Parameters.Transition.Text;
    zoneshape{n} = regs{1,n}.Shape.Text;
%     zonecentre(n,1) = str2double(regs{1,n}.Center.Attributes.Lat);
%     zonecentre(n,2) = str2double(regs{1,n}.Center.Attributes.Lon);
%     if strcmp( 'CIRCLE',zoneshape{n});
%         zoneradius(n) = str2double(regs{1,n}.CircleRadius.Text);
%     else
%         zonepoints{n} = str2num(regs{1,n}.MapPoints.Text); %#ok<ST2NM>
%         Y{n} = zonepoints{n}(:,1) ;
%         X{n} = zonepoints{n}(:,2) ;
%     end
% from kml file - shape lats/lons
    x{n} = s(n).X;
    y{n} = s(n).Y;
end;

% Check the xml and KML files are internally consistent
% find LRM etc
isLRM = strncmp(zonetrans,'TO_LRM__',6);
isSAR = strncmp(zonetrans,'TO_SAR__',6);
isSIN = strncmp(zonetrans,'TO_SIN__',6);
isCAL = strncmp(zonetrans,'TL_ALLCS',6);
% find shape type
isPOLYGON = strcmp(zoneshape,'POLYGON');
isCIRCLE = strcmp(zoneshape,'CIRCLE');
isRECTANGLE = strcmp(zoneshape,'RECTANGLE');

% find where shape and mode match (rectangle = special case of polygon)
klrm = find(isLRM & (isPOLYGON | isRECTANGLE)) ;
ksar = find(isSAR & (isPOLYGON | isRECTANGLE));
ksin = find(isSIN & (isPOLYGON | isRECTANGLE));
kcal = find(isCAL & (isPOLYGON | isRECTANGLE));

klrm_circ = find(isLRM & isCIRCLE) ;
ksar_circ = find(isSAR & isCIRCLE );
ksin_circ = find(isSIN & isCIRCLE );
kcal_circ = find(isCAL & isCIRCLE );

% checks have everything balanced
if (sum(isSAR) + sum(isLRM) + sum(isSIN)  + sum(isCAL)) ~= numregs, error('types do not sum'), end;
if (sum(isPOLYGON) + sum(isCIRCLE) + sum(isRECTANGLE)) ~= numregs, error('shapes do not sum'), end;
if (sum(length(klrm)) + sum(length(ksar)) + sum(length(ksin)) + sum(length(kcal)) + ...
        sum(length(klrm_circ)) + sum(length(ksar_circ)) + ...
        sum(length(ksin_circ)) + sum(length(kcal_circ))) ~= numregs;
    error('shapes/types do not sum');
end;

% Reset names of transitions to simplify
zonetrans(isSAR) = {'SAR'};
zonetrans(isSIN) = {'SIN'};
zonetrans(isLRM) = {'LRM'};
zonetrans(isCAL) = {'CAL'};

% We save zoneid (name), zonetrans (mode) and zone priority
Z(:,1) = zoneid;
Z(:,2) = zonetrans;
Z(:,3) = zoneprioritycell;

% cell2struct converts cell array to nreg x 1 structure array
Z = cell2struct( Z,{'name','mode','priority'},2);
save(['zonesNameMode_' versionmask '.mat'], 'Z');

clear regs zoneid zonetrans zoneprioritycell zoneshape
clear isSAR isLRM isSIN isCAL isPOLYGON isCIRCLE isRECTANGLE
clear klrm ksar ksin kcal klrm_circ ksar_circ ksin_circ kcal_circ
%% Find & replace the Arctic and Antarctica polygons - these are the same in all mask versions
load('seasMask_polar.mat') % loads polar masks

% find CYSSAR01 and CYSSAR02 (Arctic and Antarctic) polygons
Names = {Z(:).name};
knorth = find(strncmp(Names,'CYSSAR01',8));
ksouth = find(strncmp(Names,'CYSSAR02',8));
narct = length(knorth);

% If we only have a single version - nothing else to do
if narct==1
  disp('Only a single Polar mask value included, replacing with saved polar masks');
else
  if length(ksouth)~=narct || narct~=24
    error(['There are ' num2str(narct) ' entries for CYSSAR01 and '...
      num2str(length(ksouth)) ' entries for CYSSAR02: '...
      'There should be 24 for each']);
  end
  
  % remove all but first of the North & South polar records from s, Z, x & y
  Z([knorth(2:end) ksouth(2:end)]) = [];
  s([knorth(2:end) ksouth(2:end)]) = [];
  x([knorth(2:end) ksouth(2:end)]) = [];
  y([knorth(2:end) ksouth(2:end)]) = [];
  
  % find where Arctic and Antartctic Polygons are in shortened structures.
  knorth = knorth(1); ksouth = ksouth(1);
  % Assume contiguous
  if ksouth>knorth, ksouth=ksouth-23; else knorth=knorth-23; end
  
  % reset numregs
  numregs = numregs-46;
end
%% split polygon 'CP40_003-00' into two
ksplit=find(strncmp({s(:).name},'CP40_003-00',11));

% ONLY if it exists (does not in the earlier versions of the mask)
if ~isempty(ksplit) ;
  x1 = x{ksplit}; y1 = y{ksplit}; % create copy of shape to split
  x1(y1 > 80) = []; y1(y1 > 80) = []; % remove values in Arctic
  x2 = x1(x1 > 0); y2 = y1(x1 > 0); % create copy of values E of Greenwich
  y1 = y1(x1 < 0); x1 = x1(x1 < 0); % and copy with values W of Greenwich
  x1(end+1) = x1(1); y1(end+1) = y1(1); % close shapes
  x{ksplit} = x1; y{ksplit} = y1; % Replace original shape with 'W' shape
  x{numregs+1} = x2; y{numregs+1} = y2; % Add 'E' shape to end of shapes
  % Add a new entry to structures 's' and 'Z' for this 'new' region
  Z(end+1).name = {'CP40_003'}; % Adds new entry to end of Z
  Z(end).mode = {'SAR'}; % New entry is now end!
  Z(end).priority = 4;
  s(end+1).Geometry ={'Polygon'}; % Adds new entry to end of s
  s(end).name = 'CP40_003-02';
  s(ksplit).name = 'CP40_003-01'; % Rename original kml segment correctly
end

% reset numregs
numregs = length(s);

%% Build time dependant mode mask with distinction of modes

% Mode and priority are the same for all times
Mode = {Z(:).mode};
p = [Z(:).priority];
% Take copies of the polygon coordinates
X = x;
Y = y;
% For each of the polar defined segments
%  - there are 24 (bi-monthly) in each of Arctic and Antarctic
for k=1:24
  
  pole = polar{k}; %#ok<USENS> % Find coordinates for mask for this time
  % Split the coordinates at NaNs to give North and South
  [yp,xp] = polysplit(pole(2,:),pole(1,:));
  % Replace the remainig CYSSAR01 record with North polar region coordinates
  X{knorth} = xp{1}'; Y{knorth} = yp{1}';
  % Replace the remaining CYSSAR02 record with South polar region coordinates
  X{ksouth} = xp{2}'; Y{ksouth} = yp{2}';
  
  % Save the set of Mask areas for this time to structure array
  seasMaskModes.(['seas_' sprintf('%1.2d',k)]).X = X;
  seasMaskModes.(['seas_' sprintf('%1.2d',k)]).Y = Y;
  seasMaskModes.(['seas_' sprintf('%1.2d',k)]).mode = Mode;
  seasMaskModes.(['seas_' sprintf('%1.2d',k)]).priority = p;
end
save(['seasMaskAllModes_' versionmask '.mat'],'seasMaskModes');
% -------------------------------------------------------------------------
%find LRM polygons
isLRM = strncmp(Mode,'LRM',3);
klrm = find(isLRM);
numlrm = length(klrm);
xlrm = cell(numlrm,1); ylrm = xlrm;

for i=1:numlrm
  xlrm{i} = s(klrm(i)).X;
  ylrm{i} = s(klrm(i)).Y;
end

% reset values for knorth and ksouth to account for removed records
for i=1:numlrm
  if klrm(i) < knorth, knorth = knorth-1; end
  if klrm(i) < ksouth, ksouth = ksouth-1; end
end
%% build seasonal combined mask for all seasons
for k=1:24
  % Start with the seasonal masks - Polar regions replaced
  X = seasMaskModes.(['seas_' sprintf('%1.2d',k)]).X;
  Y = seasMaskModes.(['seas_' sprintf('%1.2d',k)]).Y;
  % remove the LRM masks
  X(klrm) = [];
  Y(klrm) = [];
  
  % compute unions between all non-LRM polygons
  K = cell(numregs-numlrm,1);
  for i=1:numregs-numlrm
    for j=1:numregs-numlrm
      
      if strcmp(versionmask, '3_1') && i==1 && j ==71; % problem here as crosses dateline
        % ----- HMS - will need to check if this is still 71!
        X{j}(X{j} < 0) = X{j}(X{j} < 0) + 360 ;
      end
      % Find out of these polygons intersect
      [xi,~] = polybool('&',X{i},Y{i},X{j},Y{j});
      % If there is an intersection - we include in the union later
      if ~isempty(xi), K{i} = [K{i} j]; end
    end
  end
  
  xi = [];
  yi = [];
  for i=1:numregs-numlrm
    for j=1:length(K{i}) % For each pair of intersecting regions
      % Generate the union
      [xi,yi] = polybool('union',xi,yi,X(K{i}(j)),Y(K{i}(j)));
    end
  end
  % Convert the polygon contours back to clockwise vertices
  [xi,yi] = poly2cw(xi,yi);
  
  % remove continental ice sheets areas where LRM is operated
  for i=1:length(xi)
    for j=1:length(xlrm)
      [xi{i},yi{i}] = polybool('minus',xi{i},yi{i},xlrm{j},ylrm{j});
    end
  end
  % Convert cell arrays back to column vectors
  [yi,xi] = polyjoin(yi,xi);
  xi(end+1) = NaN; %#ok we need to add NaN so the isNaN function performs as expected
  yi(end+1) = NaN; %#ok
  % Save single 2D array of lats & lons to structure
  seasMask.(['seas_' sprintf('%1.2d',k)]) = [xi';yi'];
end

save(['seasMask_noLRM_' versionmask '.mat'],'seasMask');

% Convert polygons to faces matrix (this will just show last 'season'
[f,v] = poly2fv(xi,yi);
% Display the patch.
patch('Faces', f, 'Vertices', v, 'FaceColor', 'blue', ...
  'EdgeColor', 'none');
axis off, axis equal

