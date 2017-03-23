function [C,ymedian] = hist2d(z,xylim,nx,ny)
%
% xylim = [xmin xmax ymin ymax]
% nx is the number of bins in x
% ny is the number of bins in y

[~,n2] = size(z);
if n2 ~= 2
    error('hist2d:wrongDimension','Size of input matrix must be [n,2]')
end

if length(xylim) ~= 4
    error('hist2d:wrongRangeFormat','Range must be a vector of length 4')
end

if(xylim(1) > xylim(2) || xylim(3) > xylim(4))
    error('hist2d:wrongRangeOrder','Range order must be [xmin xmax ymin ymax]')
end

dx = (xylim(2)-xylim(1))/nx;
dy = (xylim(4)-xylim(3))/ny;
edges = {xylim(1):dx:xylim(2),xylim(3):dy:xylim(4)};
if ~all(size(z))
   z = NaN(100,2);
end
[N,C] = hist3(z,'Edges',edges);

N(N == 0) = NaN;
%pcolor(C{1},C{2},N'),shading flat
%contourf(C{1},C{2},N',20,'LineColor','none');
N(isnan(N)) = min(N(:))-(max(N(:))-min(N(:)))/5;
imagesc(flipud(C{1}),flipud(C{2}),N');
set(gca,'YDir','normal')
%set(hi,'AlphaData',~isnan(N')) % out with segmentation fault
if isnan(max(N(:)))
    caxis([0 10])
else
    caxis([0 max(N(:))])
end

% calculate the median of the SSH anomaly for each SWH bin
ymedian = NaN(1,length(edges{1}));
for i=1:length(edges{1})-1
    K = z(:,1) >= edges{1}(i) & z(:,1) < edges{1}(i+1) & ...
        z(:,2) >= edges{2}(1) & z(:,2) < edges{2}(end);
    ymedian(i) = nanmedian(z(K,2));
end
        
colormap([1 1 1;linspecer]); % bins with 0 counts will be drawn as white



