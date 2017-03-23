clear all
% reads ground tracks as downloaded from ESA ftp site
% From https://earth.esa.int/web/guest/-/ground-tracks-7209:
% FTP server: ftp://calval-pds.cryosat.esa.int
% Username: ground
% Password: tracks


cd /noc/users/cryo/QCV_Cryo2/code/groundTracks/2016_2017


fid = fopen('groundtrack_20160701T120000_20170704T170000_0001.EEF_2sec.txt','r');
fgets(fid); fgets(fid); fgets(fid);
A = fscanf(fid,'%u- %u- %u_ %u: %u: %u %*c %u %f %*c %f %f',[1,Inf]);
A = reshape(A,10,length(A)/10);
A(7,:) = [];
save groundtrack_20160701T120000_20170704T170000_0001.EEF_2sec.mat A -mat
fclose(fid);