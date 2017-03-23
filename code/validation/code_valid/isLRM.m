function k = isLRM(x,y,mask)
% Points over LRM mode areas using the ESA mode mask
%
% USAGE : k = isLRM(x,y,gshhs_x)
%
% INPUT :    x : vector containing the longitude of points to be tested.
%            y : vector containing the latitude of points to be tested
%            mask : 2 x N matrix storing the mode mask in the form of
%                   closed and non-intersecting polygons [xv;yv].
%                   The polygons represent non-LRM mode areas
% 
% OUTPUT :   k : vector the same size as x and y containing 0's and 1's. A
%                value of 1 indicates that a point is over an LRM area
%
%
% Author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
%

x = x(:)';
y = y(:)';

if(max(x) > 180)
    error('isLRM:chkCoordConvention',['the coordinate conventions for x' ...
        'and xh are not consistent (-180 to 180 and 0 to 360)'])
end


%No polygon dateline-splitting has been applied, hence some polygons span
%longitudes up to 195. We must account for this
k = ~inPolygonMex([x;y],mask);
