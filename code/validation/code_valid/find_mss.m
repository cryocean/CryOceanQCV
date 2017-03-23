function mss = find_mss(xtg,ytg)
% path1 = '/scratch/general/cryosat/validation_data_2/';
% UPDATED JAN 2017
path1 = '/noc/mpoc/cryo/cryosat/validation_data_2/';

fid = netcdf.open([path1 'DTU10MSS_1min.nc'],'NC_NOWRITE');
x = netcdf.getVar(fid,0,'double');
y = netcdf.getVar(fid,1,'double');
if xtg < 0
    xtg = xtg + 360;
end

ybottom = ytg - 5;
ytop = ytg + 5;
ky = find(y > ybottom,1,'first');
ny = find(y < ytop,1,'last') - ky + 1;
if xtg < 355 && xtg > 5
    xleft = xtg - 5;
    xright = xtg + 5;
    kx = find(x > xleft,1,'first');
    nx = find(x < xright,1,'last') - kx + 1;
elseif xtg >= 355
    xleft1 = xtg -5;
    xright1 = 360;
    xleft2 = 0;
    xright2 = mod(xtg+5,360);
    kx = find(x > xleft1,1,'first');
    nx = find(x < xright1,1,'last') - kx + 1;
    kx2 = find(x >= xleft2,1,'first');
    nx2 = find(x < xright2,1,'last') - kx2 + 1;
    x2 = netcdf.getVar(fid,0,kx2-1,nx2,'double')+360;
elseif xtg <= 5
    xleft1 = 0;
    xright1 = xtg + 5;
    xleft2 = mod(xtg-5,360);
    xright2 = 360;
    kx = find(x >= xleft1,1,'first');
    nx = find(x < xright1,1,'last') - kx + 1;
    kx2 = find(x > xleft2,1,'first');
    nx2 = find(x < xright2,1,'last') - kx2 + 1;
    x2 = netcdf.getVar(fid,0,kx2-1,nx2,'double')-360;
end

z = netcdf.getVar(fid,2,[kx-1,ky-1],[nx,ny],'double');
x = netcdf.getVar(fid,0,kx-1,nx,'double');
y = netcdf.getVar(fid,1,ky-1,ny,'double');
Z = z;
if exist('nx2','var')
    z2 = netcdf.getVar(fid,2,[kx2-1,ky-1],[nx2,ny],'double');
    Z = zeros(length(x)+length(x2),length(y));
    Z(1:length(x),:) = z;
    Z(length(x)+1:end,:) = z2;
    x = [x;x2];
    [x,ix] = sort(x);
    Z = Z(ix,:);
end

Z = Z.*netcdf.getAtt(fid,2,'scale_factor');
netcdf.close(fid);
mss = interp2(x,y,Z',xtg,ytg);





