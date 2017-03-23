%% |fastscatterm| documentation
% The |fastscatterm| function places color-scaled point markers on map coordinates. This is 
% a much faster version of the Mapping Toolbox's |scatterm| function, adapted
% from <http://www.mathworks.com/matlabcentral/fileexchange/47205 Aslak Grinsted's |fastscatter|>. 
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
% |fastscatterm(lat,lon,C)| places markers at geocoordinates |lat|, |lon|,
% color-scaled to values in array |C|. 
%
% |fastscatterm(...,MarkerType)| specifies a |MarkerType| as |'+'|, |'x'|, 
% |'o'|, etc. Default |MarkerType| is |'.'|.  
%
% |fastscatterm(...,'MarkerProperty',MarkerValue)| specifies any MarkerSpec 
% preferences as propery name-value pairs.  For example,
% |fastscatterm(lat,lon,z,'markersize',30)|. 
%
% |h = fastscatterm(...)| returns a handle |h| of plotted mesh object. 
% 
%% Requirements 
% This function requires Matlab's Mapping Toolbox. 

%% Example 1: Comparing |scatterm| and |fastscatterm| 
% Here we compare Matlab's inbuilt |scatterm| with |fastscatterm|.  First
% create 200,000 random data points: 

N = 200000; 
lat = 20*randn(N,1); 
lon = 15 + 10*randn(N,1); 
z = 8+7*cosd(lat)+2*sind(lon)+randn(N,1); 

%% 
% Plot with Matlab's |scatterm| function and use |tic| |toc| to calculate
% plotting time: 

figure
worldmap('africa') 

tic 
scatterm(lat,lon,15,z,'filled') 
scattermtime = toc 

caxis([7 20])

%% 
% Make an equivalent plot with |fastscatterm|: 

figure
worldmap('africa') 

tic 
fastscatterm(lat,lon,z) 
fastscattermtime = toc 

caxis([7 20])

%% 
% The difference in plotting time for |scatterm| versus |fastscatterm| is
% quite staggering. See: 

scattermtime/fastscattermtime 

%% 
% Matlab's inbuilt |scatterm| function requires more than 40 times more processing
% time than |fastscatterm| for 200,000 points.  

%% Example 2: Formatting: 
% Create 10 points of example data: 

N = 10; 
lat = 20*randn(N,1); 
lon = 15 + 10*randn(N,1); 
z = 8+7*cosd(lat)+2*sind(lon)+randn(N,1); 

%% 
% Plot the example data as big fat plus signs: 

figure
worldmap('africa') 
fastscatterm(lat,lon,z,'+','markersize',30,'linewidth',15)

%% Author Info
% This function was mostly written by <http://www.glaciology.net/ Aslak Grinsted> in 2014 based on an idea by
% Boris Babic (http://www.mathworks.com/matlabcentral/newsreader/view_thread/22966). 
% In October 2015, <http://www.chadagreene.com Chad A. Greene> of the University of Texas at Austin's Institute for 
% Geophysics merely added some coordinate transformation bits, some error checks, 
% and a little bit of documentation. 