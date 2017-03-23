function [x,y,z,t] = read_tgData_ascii(nstation)

%UPDATED JAN 2017
path1 = '/noc/mpoc/cryo/cryosat/validation_data_2/tg_data_hourly/';

% path1 = '/scratch/general/cryosat/validation_data_2/tg_data_hourly/';
nstation = sprintf('%03d',nstation);
fid = fopen([path1,'h',nstation,'.dat'],'r');
fmt = repmat('%d ',1,12);
h1 = 0:11;
h2 = 12:23;

fy = 1;
fx = 1;
head = fgetl(fid);
y = str2num(head(22:23));
yd = str2num(head(25:28))/60;
if strcmp(head(29),'S')
    fy = -1;
end
y = fy*(y + yd);

x = str2num(head(37:39));
xd = str2num(head(41:44))/60;
if strcmp(head(45),'W')
    fx = -1;
end
x = fx*(x + xd);
fclose(fid);


fid = fopen([path1,'h',nstation,'.dat'],'r');
z = [];
t = [];
while ~feof(fid)
    head = fgetl(fid); %#ok<*NASGU>
    s = textscan(fid,['%20c ',fmt]);
    
    ti = s{1};
    ti = ti(:,12:end);
    if ~feof(fid)
        ti = ti(1:end-1,:);
    end
    
    Y = str2num(ti(:,1:4)); %#ok<*ST2NM>
    M = str2num(ti(:,5:6));
    D = str2num(ti(:,7:8));
    H = str2num(ti(:,9));
    Y = Y(:,ones(12,1));
    M = M(:,ones(12,1));
    D = D(:,ones(12,1));
    H = H(:,ones(12,1));
    H1 = h1(ones(size(H,1),1),:);
    H2 = h2(ones(size(H,1),1),:);
    K = H == 1;
    H(K) = H1(K);
    H(~K) = H2(~K);
    ti = datenum(Y,M,D,H,0,0);
    s(1) = [];
    s = cellfun(@double,s,'unif',0);
    s = cell2mat(s);
    
    
    ti = ti';
    s = s';
    t = [t;ti(:)]; %#ok<*AGROW>
    z = [z;s(:)];
end
fclose (fid);
z(z == 9999 | abs(z) > 1e+5) = NaN;

