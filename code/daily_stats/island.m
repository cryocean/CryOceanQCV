function k = island(x,y,gshhs_x,lake_flag)
% Points over land using GSHHS shoreline data
%
% USAGE : k = island(x,y,gshhs_x)
%
% INPUT :    x : vector containing the longitude of points to be tested.
%            y : vector containing the latitude of points to be tested
%            gsshs_x : structure containing the fields: lon,lat, and level.
%                      The structure refers to the GSHHS shoreline data
%            lake_flag : set lake_flag equal to 1 if lakes and enclosed
%                        seas are to be excluded of what we consider "land"
%                        default is 0.
% 
% OUTPUT :   k : vector the same size as x and y containing 0's and 1's. A
%                value of 1 indicates that a point is over land
%
% Note: the script uses the boundary between Antarctica grounding-line and
% ocean rather than that between Antarctica ice and ocean. This choice has
% been hard-coded but can ben modified by changing the value of level.
%
% Author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
%

if nargin < 4
    lake_flag = 0;
end

x = x(:)';
y = y(:)';
xh = gshhs_x.lon;
yh = gshhs_x.lat;
lev = gshhs_x.level;

%boundary between Antarctica grounding-line and ocean (level 6) rather than
%Antarctica ice and ocean (level 5)
if lake_flag == 0
    kl = (lev == 1 | lev == 6);
else
    kl = (lev == 1 | lev == 2 | lev == 6); %include lakes and enclosed seas
end
xh = xh(kl);
yh = yh(kl);
xh = cell2mat(xh);
yh = cell2mat(yh);
xh = xh(:)';
yh = yh(:)';

if(min(x) < 0 && min(xh) > 0)
    error('isLand:chkCoordConvention',['the coordinate conventions for x' ...
        'and xh are not consistent (-180 to 180 and 0 to 360)'])
end

%Since the Antarctic polygon spans all longitudes it is not closed, and
%thus we need to close it.
ka = find(yh <-62 & xh < 180.5);
k0 = ka(xh(ka) == min(xh(ka)));
if(length(xh) < 1e6 && length(xh)>1e4)
    % intermediate resolution
    kf = ka(xh(ka) > 179.5 & xh(ka) < 179.61 & yh(ka) > -84.33 & yh(ka) < -84.3);
elseif length(xh)>1e6
    kf = ka(xh(ka) == max(xh(ka))); %high resolution
else % coarse resolution
    kf = ka(xh(ka) > 167.14 & xh(ka) < 167.16 & yh(ka) > -83.25 & yh(ka) < -83.22);
end
xh = [xh(1:k0-1) -180 xh(k0:kf) 180 -180 xh(kf+1:end)];
yh = [yh(1:k0-1) -90 yh(k0:kf) -90 -90 yh(kf+1:end)];

%No polygon dateline-splitting has been applied, hence some polygons span
%longitudes up to 195. We must account for this
kr = x <= -165;
k = inPolygonMex([x;y],[xh;yh]);
if sum(kr) ~= 0
    kp = inPolygonMex([x(kr)+360;y(kr)],[xh;yh]);
    k(kr) = or(k(kr),kp);
end







