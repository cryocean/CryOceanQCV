function h = fastscatterm_NOC(lat,lon,C,varargin)
% fastscatterm places color-scaled point markers on map coordinates. This is 
% a much faster version of the Mapping Toolbox's scatterm function, adapted
% from Aslak Grinsted's fastscatter 
% (http://www.mathworks.com/matlabcentral/fileexchange/47205). 
% modified by Chris Banks, 19 June 2017 so will work with TDS data
% 
%% Syntax
% 
%  fastscatterm(lat,lon,C)
%  fastscatterm(...,MarkerType)
%  fastscatterm(...,'MarkerProperty',MarkerValue) 
%  h = fastscatterm(...)
% 
%% Description 
% 
% fastscatterm(lat,lon,C) places markers at geocoordinates lat, lon,
% color-scaled to values in array C. 
%
% fastscatterm(...,MarkerType) specifies a MarkerType as '+', 'x', 
% 'o', etc. Default MarkerType is '.'.  
%
% fastscatterm(...,'MarkerProperty',MarkerValue) specifies any MarkerSpec 
% preferences as propery name-value pairs.  For example,
% fastscatterm(lat,lon,z,'markersize',30). 
%
% h = fastscatterm(...) returns a handle h of plotted mesh object. 
% 
%% Requirements 
% This function requires Matlab's Mapping Toolbox. 
% 
%% Examples
% % With this sample data of 200,000 random points: 
% 
% N = 200000; 
% lat = 20*randn(N,1); 
% lon = 15 + 10*randn(N,1); 
% z = 8+7*cosd(lat)+2*sind(lon)+randn(N,1); 
% 
% % Plot with fastscatterm: 
% fastscatterm(lat,lon,z) 
% 
% % Alternatively, plot as big fat plus signs: 
% fastscatterm(lat,lon,z,'+','markersize',30,'linewidth',15)
% 
%% Author Info
% This function was mostly written by Aslak Grinsted in 2014 based on an idea by
% Boris Babic (http://www.mathworks.com/matlabcentral/newsreader/view_thread/22966). 
% In October 2015, Chad Greene of the University of Texas at Austin's Institute for 
% Geophysics merely added some coordinate transformation bits, some error checks, 
% and a little bit of documentation. 
% 
% See also scatterm, scatter, plotm, and fastscatter. 

%% Error checks:

% % % % % assert(license('test','map_toolbox')==1,'The fastscatterm function requires Matlab''s Mapping Toolbox.') 
% % % % % assert(ismap(gca)==1,'You must first initialize a map before calling fastscatterm.') 
assert(isequal(size(lat),size(lon))==1,'Inputs lat and lon must be the same size.') 
assert(length(C)==length(lat),'Color data C must have dimensions corresponding to the dimensions of lat and lon.')

%% Set defaults: 

marker = '.'; 
markersize = 2; 

%% Parse Inputs: 

if nargin>3
    if any(strcmp(varargin{1},{'.','o','x','+','*','s','d','v','<','^','<','p','h'}))
    	marker=varargin{1};
        varargin(1)=[];
    end
    
    tmp = strcmpi(varargin,'markersize'); 
    if any(tmp)
        markersize = varargin{find(tmp)+1}; 
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
    end
end

%% Coordinate transformations: 

[X,Y] = mfwdtran(lat,lon); 

% %h/t to Boris Babic for the method. see http://www.mathworks.com/matlabcentral/newsreader/view_thread/22966
% h=mesh([X(:) X(:)]',[Y(:) Y(:)]',zeros(2,numel(X)),'mesh','column','marker',marker,'cdata',[C(:) C(:)]',varargin{:});
ix=find(~isnan(C+X+Y));
if mod(length(ix),2)==1
    ix(end+1)=ix(end);
end
ix=reshape(ix,2,[]);

h=mesh(X(ix),Y(ix),zeros(size(ix)),'marker',marker,'cdata',C(ix),'markersize',markersize,...
    'edgecolor','none','markeredgecolor','flat','facecolor','none',varargin{:});

%% Set view: 

view(2)
grid off

%% Clean up: 

if nargout==0
    clear h
end