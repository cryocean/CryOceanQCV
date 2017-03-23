clear all
fn = dir('/noc/users/cryo/QCV_Cryo2/code/mode_mask/*.shp');
s = shaperead(fn(1).name,'UseGeoCoords',true);
[y,x] = polysplit([s(:).Lat],[s(:).Lon]);

% find LRM polygons
xlrm = cell(length(x),1);
ylrm = cell(length(x),1);
ct = 0;
klrm = zeros(length(x),1);
for i=1:length(x)
    if(strfind(s(i).Transition,'LRM'))
        ct = ct+1;
        klrm(ct) = i;
        xlrm{ct} = s(i).Lon;
        ylrm{ct} = s(i).Lat;
    end
end
klrm(klrm == 0) = [];
k = cellfun(@isempty,xlrm);
xlrm(k) = [];
ylrm(k) = [];
x(klrm) = [];
y(klrm) = [];

% we are interested in LRM and !LRM
xlrm1 = xlrm{1}; % LRM polygon in Antarctica
ylrm1 = ylrm{1};
% extend LRM polygon down to -90 and close it
k = find(ylrm1 < -84); % there are only two points with latitude < -84
xlrm{1} = [xlrm1(1:k(1)) 180 180 -180 -180 xlrm1(k(2):end)];
ylrm{1} = [ylrm1(1:k(1)) ylrm1(k(1)) -90 -90 ylrm1(k(2)) ylrm1(k(2):end)];

% extend no-LRM Antarctica polygon down to -90 and close it
x1 = x{1};
y1 = y{1};
kr = find(x1 == max(x1));
kl = find(x1 == min(x1));
x1 = [x1(1:kr) 180 180 -180 -180 x1(kl:end)];
y1 = [y1(1:kr) y1(kr) -90 -90 y1(kl) y1(kl:end)];
x{1} = x1;
y{1} = y1;

% extend no-LRM Arctic polygon up to 90 and close it
x1 = x{18};
y1 = y{18};
kl = find(x1 == min(x1));
kr = kl-1;
x1 = [x1(1:kr) 180 180 -180 -180 x1(kl:end)];
y1 = [y1(1:kr) y1(kr) 90 90 y1(kl) y1(kl:end)];
[x1,y1] = poly2cw(x1,y1);
x{18} = x1;
y{18} = y1;

% split polygon 71 into two
x1 = [x{71}(1) x{71}(4) 180 180 x{71}(1)];
y1 = [y{71}(1) y{71}(4) y{71}(4) y{71}(1) y{71}(1)];
x2 = [x{71}(2) -180 -180 x{71}(3) x{71}(3)];
y2 = [y{71}(2) y{71}(2) y{71}(3) y{71}(3) y{71}(2)];
x{71} = x1;
y{71} = y1;
x{length(x)+1} = x2;
y{length(y)+1} = y2;

% remove polygons 2 and 19 (they are duplicated of 1 and 18 !?)
x([2 19]) = [];
y([2 19]) = [];

% compute unions between all no-LRM polygons
K = cell(length(x),1);
for i=1:length(x)
    for j=1:length(x)
        [xi,~] = polybool('&',x{i},y{i},x{j},y{j});
        if ~isempty(xi)
            K{i} = [K{i} j];
        end
    end
end

xi = [];
yi = [];
for i=1:length(x)
    for j=1:length(K{i})
        [xi,yi] = polybool('union',xi,yi,x(K{i}(j)),y(K{i}(j)));
    end
end
[xi,yi] = poly2cw(xi,yi);

% remove continental ice sheets areas where LRM is operated
for i=1:length(xi)
    [xi{i},yi{i}] = polybool('minus',xi{i},yi{i},xlrm{1},ylrm{1});
    [xi{i},yi{i}] = polybool('minus',xi{i},yi{i},xlrm{2},ylrm{2});
end
[yi,xi] = polyjoin(yi,xi);

% [f,v] = poly2fv(xi,yi);
% % Display the patch.
% patch('Faces', f, 'Vertices', v, 'FaceColor', 'blue', ...
%  'EdgeColor', 'none');
% axis off, axis equal

        
    

