
function dyn_argo_profiles

% ---------------------------- set paths ----------------------------------
% UPDATED JAN 2017
path1 = '/noc/mpoc/cryo/cryosat/validation_data_2/EN4_TS_profiles/';

% path1 = '/scratch/general/cryosat/validation_data_2/EN4_TS_profiles/';
% fn = 'EN.4.1.1.f.profiles.g10.'; % old form
fn = 'EN.4.2.0.f.profiles.g10.'; % final form

addpath(genpath('/noc/users/cryo/QCV_Cryo2/code/validation/code_valid/gsw_matlab_v3_04/'))
zref = 1000;
% -------------------------------------------------------------------------

year0 = 2014;
tnow = datenum(date,'dd-mmm-yyyy');
[Y,~,~] = datevec(tnow);
yearf = Y;

for nyear=year0:yearf;
    for nmonth=1:12;
        disp([nyear,nmonth])
        if exist([path1 'dataArgo_' num2str(nyear) '_' num2str(nmonth,'%02i'), ...
                '_' num2str(zref) '.mat'],'file')==0 && ...
            exist(sprintf([path1 fn '%i%02i.nc'],nyear,nmonth),'file') == 2
            
            % create array of empty structures with all the fieldnames. There
            % will be one structure per Argo profile and per month
            dataArgo = repmat(struct('id',[],'time',[],'lon',[],'lat',[],...
                'depth',[],'T',[],'S',[],'steric',[]),10000,1);
            
            % read data
            ncid = netcdf.open(sprintf([path1 fn '%i%02i.nc'],nyear,nmonth),...
                'NC_NOWRITE');
            
            varid = netcdf.inqVarID(ncid,'LATITUDE');
            fillVal = netcdf.getAtt(ncid,varid,'_fillvalue');
            lat = netcdf.getVar(ncid,varid,'double');
            lat(lat == fillVal) = NaN;
            
            varid = netcdf.inqVarID(ncid,'LONGITUDE');
            fillVal = netcdf.getAtt(ncid,varid,'_fillvalue');
            lon = netcdf.getVar(ncid,varid,'double');
            lon(lon == fillVal) = NaN;
            
            varid = netcdf.inqVarID(ncid,'DEPH_CORRECTED');
            fillVal = netcdf.getAtt(ncid,varid,'_fillvalue');
            depth = netcdf.getVar(ncid,varid,'double');
            depth(depth == fillVal) = NaN;
            
            varid = netcdf.inqVarID(ncid,'JULD');
            fillVal = netcdf.getAtt(ncid,varid,'_fillvalue');
            time = netcdf.getVar(ncid,varid,'double') + datenum(1950,1,1); %UTC
            time(time == fillVal) = NaN;
            
            varid = netcdf.inqVarID(ncid,'POTM_CORRECTED'); % potential temp
            fillVal = netcdf.getAtt(ncid,varid,'_fillvalue');
            T = netcdf.getVar(ncid,varid,'double');
            T(T == fillVal) = NaN;
            
            varid = netcdf.inqVarID(ncid,'PSAL_CORRECTED');
            fillVal = netcdf.getAtt(ncid,varid,'_fillvalue');
            S = netcdf.getVar(ncid,varid,'double');
            S(S == fillVal) = NaN;
            
            % read quality control flags and other info
            varid = netcdf.inqVarID(ncid,'POSITION_QC');
            position_qc = netcdf.getVar(ncid,varid); %if '4' reject entire profile
            
            varid = netcdf.inqVarID(ncid,'PROFILE_POTM_QC');
            profile_T_qc = netcdf.getVar(ncid,varid); %if '4' reject entire T profile
            
            varid = netcdf.inqVarID(ncid,'PROFILE_PSAL_QC');
            profile_S_qc = netcdf.getVar(ncid,varid); %if '4' reject entire S profile
            
            varid = netcdf.inqVarID(ncid,'POTM_CORRECTED_QC');
            potm_qc = netcdf.getVar(ncid,varid);
            
            varid = netcdf.inqVarID(ncid,'PSAL_CORRECTED_QC');
            psal_qc = netcdf.getVar(ncid,varid);
            
            varid = netcdf.inqVarID(ncid,'PROJECT_NAME');
            project = netcdf.getVar(ncid,varid);
            
            varid = netcdf.inqVarID(ncid,'PLATFORM_NUMBER');
            id = netcdf.getVar(ncid,varid);
            
            netcdf.close(ncid); % close netcdf file
            
            % find Argo floats (we want to use only Argo floats)
            project = cellstr(project');
            kargo = cellfun(@strfind,project,repmat({'ARGO'},length(project),...
                1),'unif',0);
            kargo = cellfun(@(x) ~isempty(x),kargo);
            
            
            % Reject profiles and levels based on qc flags
            ki = find(position_qc ~= '4' & profile_T_qc ~= '4' & ...
                profile_S_qc ~= '4' & kargo);
            lat = lat(ki);
            lon = lon(ki);
            time = time(ki);
            depth = depth(:,ki);
            T = T(:,ki);
            S = S(:,ki);
            potm_qc = potm_qc(:,ki);
            psal_qc = psal_qc(:,ki);
            %project = project(ki);
            id = id(:,ki);
            T(potm_qc == '4') = NaN;
            S(psal_qc == '4') = NaN;
            % ---------------------------------- end rejection
            
           
            % find maximum depth of (last) observation
            [nr,nc] = find(~isnan(T));
            [~,mm] = unique(nc,'last'); % mm is the index of the last repeated element
            nr = nr(mm); nc = nc(mm);
            pmax = sub2ind(size(T),nr,nc);
            maxT = depth(pmax);
            
            [nr,nc] = find(~isnan(S));
            [~,mm] = unique(nc,'last'); % mm is the index of the last repeated element
            nr = nr(mm); nc = nc(mm);
            pmax = sub2ind(size(S),nr,nc);
            maxS = depth(pmax);
            % ---------------------------------- end finding
            
            % reject profiles that have no data below zref m
            kdepth = maxT > zref & maxS > zref;
            lat = lat(kdepth);
            lon = lon(kdepth);
            time = time(kdepth);
            depth = depth(:,kdepth);
            T = T(:,kdepth);
            S = S(:,kdepth);
            %project = project(kdepth);
            id = id(:,kdepth);
            knegative = find(depth < 0 | depth > 5000);
            T(knegative) = NaN;
            S(knegative) = NaN;
            depth(knegative) = NaN;
            % ---------------------------------- end rejection
            
            
            %---------- compute steric height for each profile ----------------
            hst = zeros(length(time),1);
            kgood = true(length(time),1);
            for i=1:length(hst)
                try
                    zi = depth(:,i);
                    kref = find(zi >= zref & ~isnan(T(:,i)) & ~isnan(S(:,i)),...
                        1,'first');
                    %pref = gsw_p_from_z(-zref,lat(i));
                    pi = gsw_p_from_z(-zi(1:kref),lat(i));
                    
                    %compute absolute salinity from practical salinity
                    [SA, ~] = gsw_SA_from_SP(S(1:kref,i),pi,lon(i),lat(i));
                    
                    %compute conservative temperature from potential temperature
                    CT = gsw_CT_from_pt(SA,T(1:kref,i));
                    
                    %compute steric height
                    hstmp = gsw_geo_strf_steric_height(SA,CT,pi,zref);
                    
                    hst(i) = hstmp(1);
                catch %#ok
                    kgood(i) = false;
                end
            end
            
            lat = lat(kgood);
            lon = lon(kgood);
            time = time(kgood);
            depth = depth(:,kgood);
            T = T(:,kgood);
            S = S(:,kgood);
            %project = project(kgood);
            id = id(:,kgood);
            hst = hst(kgood);
            
            % -----------------------------------------------------------------
            
            
            % store data for each profile in structure
            id = id';
            id = cellstr(id);
            idu = unique(id);
            for i=1:length(idu);
                kid = cellfun(@strcmp,id,repmat(idu(i),length(id),1));
                dataArgo(i).id = idu{i};
                dataArgo(i).time = time(kid);
                dataArgo(i).lon = lon(kid);
                dataArgo(i).lat = lat(kid);
                dataArgo(i).depth = depth(:,kid);
                dataArgo(i).T = T(:,kid);
                dataArgo(i).S = S(:,kid);
                dataArgo(i).steric = hst(kid);
            end
            
            kempty = arrayfun(@(X) isempty(X.T),dataArgo);
            dataArgo = dataArgo(~kempty); %#ok remove empty structures
            save([path1 'dataArgo_' num2str(nyear) '_' num2str(nmonth,'%02i'), ...
                '_' num2str(zref) '.mat'],...
                'dataArgo','-mat')
        end
    end
end
            
        
        
