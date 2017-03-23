function h=read_mssh(lon,lat,model)
% MSSH Mean Sea Surface Height from various models
%   H=READ_MSSH(LON,LAT,MODEL) outputs the mean Sea Surface H in m over the
%     locations specified by LON,LAT.
%    LON,LAT can be vectors of equal length (useful in the along-track case)
%     or matrices defining a geographical grid, for instance from MESHGRID.
%    MODEL is a string specifying the Mean Sea Surface (MSS) model:
%       'dtu10' = DTU2010 MSS (default)
%       'cls11' = CNES/CLS2011 MSS
% Note that the code checks if the selected model is already loaded in the
% current workspace, in a variable called DTU10MSS or CLS11MSS, to speed up
% execution in case of multiple calls
%
% V2.0 Paolo Cipollini 19/06/2013

% ==== VERSION HISTORY ====
% V2.0 19/06/2013 Paolo Cipollini
%           added CNES_CLS_2011 model
%           accounted for negative longitudes
% V1.0 08/03/2013 Paolo Cipollini

if nargin==2, model='dtu10'; end

isgood=(isfinite(lon) & isfinite(lat));
h=nan(size(lon));

switch model
    case 'dtu10'
        isloaded=evalin('caller','exist(''DTU10MSS'',''var'');');
        if isloaded
            MSS=evalin('caller','DTU10MSS');
        else
            %MSS=read_nc('/noc/bodc/lso/alt/mss/dtu10/DTU10MSS_2min.nc');
            MSS = read_nc('DTU10MSS_2min.nc');
        end
        
        MSS_lat = MSS.lat.Value;
        MSS_lon = MSS.lon.Value;
        MSSH = MSS.mss.Value';
        MSS_lon(end) = [];
        MSSH(:,end) = [];
        MSS_lon(MSS_lon > 180) = MSS_lon(MSS_lon > 180) - 360;
        [MSS_lon,ix] = sort(MSS_lon);
        MSSH = MSSH(:,ix);
        
    case 'cls11'
        
        isloaded=evalin('base','exist(''CLS11MSS'',''var'');');
        if isloaded
            MSS=evalin('base','CLS11MSS');
        else
            MSS=read_nc('/noc/bodc/lso/alt/mss/cls11/mss_cnes_cls2011.nc');
        end
        MSS_lat=MSS.NbLatitudes.Value+MSS.LatLonStep.Value(1)/2;
        MSS_lon=MSS.NbLongitudes.Value+MSS.LatLonStep.Value(2)/2;
        MSSH=MSS.Grid_0001.Value;
end

% h(isgood)=...
%     interp2(MSS_lon,MSS_lat,MSSH,mod(lon(isgood),360),lat(isgood));

h(isgood)=...
    interp2(MSS_lon,MSS_lat,MSSH,lon(isgood),lat(isgood));
return
