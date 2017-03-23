function indCol = dot2color(z,ncolor,cmax,cmin)
% linearly maps values between cmin and cmax to be used in a colormap with
% ncolor colors 
% z : vector containing data values
% ncolor : number of colors of colormap
% cmax : upper color limit for caxis([cmin cmax]) (default max(z))
% cmin : lower color limit for caxis([cmin cmax]) (default min(z))
%
% indCol : contains the indices that map each value in z to a colormap with
%         ncolor colors. Hence the color associated with each value is, for
%         instance,: cmap = jet(ncolor); col = cmap(indCol,:); 

if nargin < 2
    error('dot2color:NotEnoughInputs','Requires at least 2 inputs')
elseif (nargin > 2 && isempty(cmax) && nargin < 4)
    error('dot2color:cminNotSpecified',['cmin needs to be specified', ...
        ' if cmax is empty'])
end

if nargin == 2
    cmax = max(z);
    cmin = min(z);
elseif isempty(cmax)
    cmax = max(z);
end
if (nargin == 3 || isempty(cmin))
    cmin = min(z);
end

indCol = round((ncolor-1)/(cmax-cmin).*z+(cmax-cmin*ncolor)/(cmax-cmin));
indCol(indCol > ncolor) = ncolor;
indCol(indCol < 1) = 1;
    
    