function delta_h=dh_ellips(latitude)
% DH_ELLIPS Difference in elevation for WGS84 minus TOPEX ellipsoid 
% from ftp://sidads.colorado.edu/DATASETS/icesat/tools/idl/ellipsoid/README_ellipsoid.txt
%  Quoting:
% Since the demonstrated differences in latitude between these ellipsoids is
% so small, it is possible to approximate the change in elevation between
% ellipsoids for a particular latitude using an empirically-derived formula:
% 
%   delta_h = h2 - h1 =  -((a2 - a1) * cos(phi)^2 + (b2 - b1) * sin(phi)^2
%     where
%       phi is latitude.
%       h1 and h2 are elevations for ellipsoids 1 and 2, respectively.
%       a1 and a2 are equatorial radii of ellipsoids 1 and 2, respectively.
%       b1 and b2 are polar radii of ellipsoids 1 and 2, respectively.

a1 = 6378137 ; f1 = 1/298.257223563 ; b1 = a1*(1-f1);

a2 = 6378136.3 ; f2 = 1/298.257 ; b2 = a2*(1-f2);

delta_h =   -((a2 - a1) * cosd(latitude).^2 + (b2 - b1) * sind(latitude).^2);