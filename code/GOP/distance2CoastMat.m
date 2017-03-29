function d = distance2CoastMat(x,y,xc,yc)
deg2rad = pi/180.0;
R = 6371000;
distan = @(x1,y1,x2,y2) R.*acos(sin(y1.*deg2rad).*sin(y2.*deg2rad) + ...
            cos(y1.*deg2rad).*cos(y2.*deg2rad).*cos((x2-x1).*deg2rad));
x = x(:);
y = y(:);
xc = xc(:);
yc = yc(:);
n = length(x);
d = NaN(n,1);
for i=1:n
    k = abs(yc-y(i)) <= 1 & (abs(xc-x(i)) < 15 | abs(xc-x(i)) > 345);
    if any(k)
        dtmp = distan(xc(k),yc(k),x(i),y(i));
        d(i) = min(dtmp);
    end
end