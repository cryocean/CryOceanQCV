function dx = Distance(x1,y1,x2,y2,type)
% DX = Distance(x1,y1,x2,y2) compute the non-euclidean distance between two
% points on the Earth.
% INPUT : 
%          x1 : Longitude of point 1
%          y1 : Latitude of point 1
%          x2 : Longitude of point 2
%          y2 : Latitude of point 2
%
% OUTPUT : 
%          dx : distance in 'm'

% Note : type = 0 uses a more precisse function (default) while type = 1
% uses a less precise but faster function.


if nargin < 4
    error('Distance:NotEnoughInputs', 'Needs at least 4 inputs.');
end

if (nargin == 4 || type == 0)
Y = (sind(x2-x1).*cosd(y2)).*(sind(x2-x1).*cosd(y2)) + (cosd(y1).*sind(y2) ...
    - sind(y1).*cosd(y2).*cosd(x2-x1)).*(cosd(y1).*sind(y2) - sind(y1) ... 
    .*cosd(y2).*cosd(x2-x1));
X = sind(y1).*sind(y2) + cosd(y1).*cosd(y2).*cosd(x2-x1);
dx = 6372000.*atan2(sqrt(Y),X);
elseif (type == 1)
dx = 6372000*acos(sind(y1).*sind(y2)+cosd(y1).*cosd(y2).*cosd(x2-x1));
else
    error('Distance:WrongType', 'Type can only be 0 or 1');
end