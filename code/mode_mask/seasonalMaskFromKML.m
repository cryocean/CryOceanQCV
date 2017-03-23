%clear all
fn = dir('/noc/users/cryo/QCV_Cryo2/code/mode_mask/Seasonal*.kml');
s = kml_shapefile(fn(1).name);
xnorth = cell(24,1);
ynorth = cell(24,1);
xsouth = cell(24,1);
ysouth = cell(24,1);
for i=1:24
    xnorth{i} = s(i).X;
    ynorth{i} = s(i).Y;
    xsouth{i} = s(i+24).X;
    ysouth{i} = s(i+24).Y;
end
