function h = plotTable(strRow,strCol,strData,xy,fontSz,lWidth)
% plot a table with a custom number of columns and rows. It returns a
% handle to all annotation objects. The handle can be used to modify the
% properties of each table cell. For instance, if you want to set the
% background color of the cell (2,3) to red, you would write:
% set(h(2,3),'backGroundColor','red').
%
% USAGE :     h = plotTable(strRow,strCol,strData,xy,fontSz,lWidth)
%
% INPUT :     strRow : cell array of strings containing the labels of the
%                      rows. The first element corresponds to the second
%                      row in the table (i.e., grid(2,1) ). If rows have no
%                      labels, enter an empty array (i.e., strRow = []).
%             strCol : cell array of strings containing the labels of the
%                      columns. If strRow = [], the number of labels in
%                      strCol must coincide with the length of strData.
%             strData : cell array of strings containing the data.
%             xy : location of the lower left corner of the object in
%                  centimeters (x,y).
%             fontSz : size of the characters in points.
%             lWidth : width of edges of the table.
%
% OUTPUT :    h : handle to all annotation objects (i.e., all table cells).
% 
% Example : sc = {'test','mean','std'};
%           sr = {'test1','test2'};
%           sd = {'3.5','0.5';'3.7','0.6'};
%           h = plotTable2(sr,sc,sd,[4 3],12,1.5);
%           
% Author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)

if (~isempty(strRow) && length(strCol)-1 ~= size(strData,2) || ...
        isempty(strRow) && length(strCol) ~= size(strData,2))
    error('Number of column labels inconsistent with length of data')
end
if ~isempty(strRow) && length(strRow) ~= size(strData,1)
    error('Number of row labels inconsistent with length of data')
end
if nargin < 5
    fontSz = 10;
    lWidth = 1;
end
if nargin < 6
    lWidth = 1;
end

%wx = @(x,y) x*y*0.03+y*0.0117; % cell width
wx = @(x,y) x*y*0.023+y*0.01; % cell width
%wy = @(x) 0.0583*x+0.233;
%wy = @(x) 0.0583*x-0.01; % cell height
wy = @(x) 0.0583*x+0.02; % cell height
h = zeros(length(strRow)+1,length(strCol));
bgColor = [1 1 1];


% plot first column if necessary
if ~isempty(strRow)
    strRow = strRow(:);
    strRow = flipud(strRow);
    strRow(end+1) = strCol(1);
    strCol(1) = [];
    lenCol1 = max(cellfun(@length,strRow))+0.5;
    of = 1;
    posy = 0;
    for i=1:length(strRow)
        posCol = [xy(1) xy(2)+posy wx(lenCol1,fontSz) wy(fontSz)];
        h(i,1) = annotation('textbox','BackGroundColor',bgColor, ...
            'String',strRow{i},'units','centimeters','Position',posCol, ...
            'LineWidth',lWidth,'VerticalAlignment','middle', ...
            'HorizontalAlignment','center','fontSize',fontSz);
        posy = posy + wy(fontSz);
    end
else
    lenCol1 = 0;
    of = 0;
end

% plot remaining columns
posx = wx(lenCol1,fontSz);
for j=1:length(strCol)
    lenCol = max(max(cellfun(@length,strData(:,j))),length(strCol{j}));
    tmpData = flipud({strCol{j},strData{:,j}}');
    posy = 0;
    for i=1:length(tmpData)
        posCol = [xy(1)+posx xy(2)+posy wx(lenCol,fontSz) wy(fontSz)];
        h(i,j+of) = annotation('textbox','BackGroundColor',bgColor, ...
            'String',tmpData{i},'units','centimeters','Position',posCol, ...
            'LineWidth',lWidth,'VerticalAlignment','middle', ...
            'HorizontalAlignment','center','fontSize',fontSz);
        posy = posy + wy(fontSz);
    end
    posx = posx + wx(lenCol,fontSz);
end
h = flipud(h);