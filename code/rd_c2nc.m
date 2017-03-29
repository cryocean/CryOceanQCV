function [Vars,VAtts,GAtts,Dims]=rd_c2nc(filename,varargin)

% Read Cryosat 2 Baseline C Netcdf file into MATLAB structures
% Applies scaling and offset if CF attributes applied
%  VARS = RD_NC(FILENAME) reads into structure VARS all variables
%     from the NetCDF file specified by FILENAME (which may be a full path).
%     This allows full access to each variable and its attributes in the
%     MATLAB workspace.
%  VARS = RD_NC(FILENAME,Var_1,...Var_N) reads only variables named Var_1
%     to Var_N, if they exist in the file from the NetCDF file.
%  [VARS,GATTS] = READ_NC(FILENAME) also reads the Global Attributes
%     of the NetCDF file (into structure GATTS).
%  [VARS,GATTS,DIMS] = READ_NC(FILENAME) also separately reads the
%     dimensions (into structure DIMS).
%
%  v 0.0  Helen Snaith (NOC), Mar 2017
%  based on rd_nc v 4.0  Helen Snaith (NOC), Mar 2016
%  Changed to read to variables into one structure and variable attributes
%   to a second, to match formats as read from Cryosat 2
%   Baseline B (DBL) files
%--------------------------------------------------------------------------

% Get NC_GLOBAL
% ncglob = netcdf.getConstant('NC_GLOBAL');

% Find if we are looking for specific variables
if nargin>1
  v_names = cell(nargin-1,1);
  for i=1:nargin-1
    v_names{i} = char(varargin(i));
  end
end

% open file
ncid=netcdf.open(filename,'NC_NOWRITE');

% Get info on file content
finfo = ncinfo(filename);
ndims = length(finfo.Dimensions);
nvars = length(finfo.Variables);
ngatts = length(finfo.Attributes);
%[ndims,nvars,ngatts] = netcdf.inq(ncid);

% Setup an empty structures for attributes and Variables
% att_struc=struct('Type',[],'Length',[],'Value',[]);
% var_struc=struct('Type',[],'DimIds',[],'Attr',[]);

% Get the dimensions
if nargout>3
  for i = 1:ndims
    Dims.(finfo.Dimensions(i).Name) = finfo.Dimensions(i).Length;
  end
end

% Get the global attributes.
if (nargout>2 && ngatts>0)
  for i = 1:ngatts
    GAtts.(finfo.Attributes(i).Name) = finfo.Attributes(i).Value;
  end
else
  GAtts = [];
end

% Now for the Variables
if (nvars>0)
  for i=1:nvars
    % If we want this variable
    if ~exist('v_names','var') || any(strcmp(finfo.Variables(i).Name,v_names))
      % Save the attributes
      varname = finfo.Variables(i).Name;
      VAtts.(varname) = finfo.Variables(i);
      attNames = {VAtts.(varname).Attributes(:).Name};
      % Read the Data Values
      % Most scaled to doubles - calendar gets time vars
      if ismember('scale_factor',attNames) || ismember('calendar',attNames)
        datain = ncread(filename,finfo.Variables(i).Name);
%         fprintf('%s %d %d\n',varname,size(datain));
       % The rest need to stay as written in file (variety of integers)
      else
        datain = netcdf.getVar(ncid,i-1);
      end
      
      
      % May still need to do this
      % Check for unsigned integers
      %  In these cases - datain will wrap at intmax(signed integer) and
      %  should have intmax(unsigned integer) added to negative values
%       if ismember('Unsigned',{finfo.Variables(i).Attributes(:).Name}) &&...
%                     ~strncmp(class(datain),'u',1)
%         if min(datain(:))<0
%           tmpv = cast(datain,['u' class(datain)]);
%           ifneg = find(datain<0);
%           tmpv(ifneg) = double(datain(ifneg))+double(intmax(['u' class(datain)]));
%           datain = tmpv;
%           clear tmpv ifneg;
%         else
%           datain=cast(datain,['u' class(datain)]);
%         end
%         %If the FillValue has also been set to signed value
%         if isfield(Vatts.(finfo.Variables(i).Name),'FillValue')
%           Fill = VAtts.(finfo.Variables(i).Name).FillValue;
%           if ~strncmp(class(Fill),'u',1)
%             Need to correct FillValue as well
%             if Fill<0
%               Fill = cast(double(Fill)+double(intmax(['u' class(Fill)])),['u' class(Fill)]);
%             else
%               Fill = cast(Fill,['u' class(Fill)]);
%             end
%            VAtts.(finfo.Variables(i).Name).FillValue = Fill;
%           end
%       % Set default data (value == FillValue) - only needed for unsigned
%         bad_data = find(datain==Fill);
%       end
%       
      
      % Save data to output structure
      Vars.(varname) = datain;
    end
  end
end

netcdf.close(ncid);

% remap the MCD flag field
flds = fieldnames(Vars);
if ismember('flag_mcd_20_ku',flds)
  mcd = Vars.('flag_mcd_20_ku');
  mcda = VAtts.('flag_mcd_20_ku').Attributes;
  for j=1:length(mcda)
    switch mcda(j).Name
      case 'flag_masks'
        flag_mask = mcda(j).Value;
      case 'flag_meanings'
        flag_meanings = mcda(j).Value;
    end
  end
  flag_meanings = strsplit(flag_meanings,' ');
  for j=1:length(flag_mask)
    Vars.(['flag_mcd_' flag_meanings{j}]) = bitand(mcd,flag_mask(j));
%     fprintf (['flag_mcd_' flag_meanings{j} ': %d %d\n'], ...
%       min(Vars.(['flag_mcd_' flag_meanings{j}])), max(Vars.(['flag_mcd_' flag_meanings{j}])));
%     Vars.('flag_mcd_block_degraded')      = bitget(mcd,32-1);
%     Vars.('flag_mcd_blank_block')         = bitget(mcd,32-1);
%     Vars.('flag_mcd_orbit_propag_err')    = bitget(mcd,32-3);
%     Vars.('flag_mcd_orbit_file_change')   = bitget(mcd,32-4);
%     Vars.('flag_mcd_orbit_discontinuity') = bitget(mcd,32-5);
%     Vars.('flag_mcd_echo_saturation')     = bitget(mcd,32-6);
%     Vars.('flag_mcd_other_echo_err')      = bitget(mcd,32-7);
%     Vars.('flag_mcd_cal1_corr_miss')      = bitget(mcd,32-12);
%     Vars.('flag_mcd_cal1_corr_ipf')       = bitget(mcd,32-13);
%     Vars.('flag_mcd_doris_uso_corr')      = bitget(mcd,32-14);
%     Vars.('flag_mcd_trk_echo_err')        = bitget(mcd,32-16);
%     Vars.('flag_mcd_rx1_echo_err')        = bitget(mcd,32-17);
%     Vars.('flag_mcd_rx2_echo_err')        = bitget(mcd,32-18);
%     Vars.('flag_mcd_cal2_corr_miss')      = bitget(mcd,32-25);
%     Vars.('flag_mcd_cal2_corr_ipf')       = bitget(mcd,32-26);
%     Vars.('flag_mcd_power_scaling_err')   = bitget(mcd,32-27);
%     Vars.('flag_mcd_processing_type')     = ...
%       bitget(mcd,32-28) + 2.*bitget(mcd,32-29);
  end
end


if ~exist('Vars','var'), Vars = []; end

end