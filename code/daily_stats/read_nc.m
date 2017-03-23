function [Vars,GAtts,Dims]=read_nc(filename,varargin)

% Read generic Netcdf file into MATLAB structures
% Applies scaling and offset if CF attributes applied
%  VARS = READ_NC(FILENAME) reads into structure VARS all variables
%     from the NetCDF file specified by FILENAME (which may be a full path).
%     This allows full access to each variable and its attributes in the
%     MATLAB workspace.
%  [VARS,GATTS] = READ_NC(FILENAME) also reads the Global Attributes
%     of the NetCDF file (into structure GATTS). 
%  [VARS,GATTS,DIMS] = READ_NC(FILENAME) also separately reads the
%     dimensions (into structure DIMS).
% An optional list of variable names can be given, in which case VARS will
%   only contain the variables named in the list, if they exist in the file
%
%  v 1.0  Helen Snaith (NOC), April 2010. rev#1 P. Cipollini 23/04/2010
%  v 2.0  Helen Snaith (NOC), Sep 2012
%  v 3.0  Helen Snaith (NOC), Mar 2014
%--------------------------------------------------------------------------
% This is a contribution to the COASTALT Project 
% ESA/ESRIN Contract No. 21201/08/I-LG, http://www.coastalt.eu
%--------------------------------------------------------------------------

% VERSION HISTORY
%  v. 1.0  Helen Snaith (NOC), April 2010
%  checked (with some variable renaming) and made into a function by cipo
%  v. 2.0 Helen Snaith (NOC), Sep 2012
%  more generic variable and attribute renaming to ensure matlab-safe field
%  names and deals with character variables
%  v. 3.0 Helen Snaith (NOC), Mar 2014
%  extend variable re-naming to dimension names
%  add capability to only extract limited variable set

% Get NC_GLOBAL
ncglob = netcdf.getConstant('NC_GLOBAL');

% Find if we are looking for speciifc variables
if nargin>1
  for i=1:nargin-1
    v_names{i} = char(varargin(i));
  end
end

% open file
ncid=netcdf.open(filename,'NC_NOWRITE');

% Get info on file content
[ndims,nvars,ngatts] = netcdf.inq(ncid);

% Setup an empty structures for attributes and Variables
att_struc=struct('Type',[],'Length',[],'Value',[]);
var_struc=struct('Type',[],'DimIds',[],'Attr',[],'Value',[]);

% Get the dimensions
for dimid = 1:ndims
  [dimName, dimLength] = netcdf.inqDim(ncid,dimid-1);
  dimName = matlab_sanitize_name(dimName);
  Dims.(dimName).Length = dimLength;
end

% Get the global attributes.
if (ngatts>0)
  for attnum = 1:ngatts
    attName = netcdf.inqAttName(ncid,ncglob,attnum-1);
    [att_struc.Type,att_struc.Length] = netcdf.inqAtt(ncid,ncglob,attName);
    att_struc.Value = netcdf.getAtt(ncid,ncglob,attName);
    attName=matlab_sanitize_name(attName);
    GAtts.(attName) = att_struc;
  end
else
  GAtts = [];
end

% Now for the Variables
if (nvars>0)
  for varid=1:nvars
    Fill = NaN;
    Missing = NaN;
    scale = 1;
    offset=0;
    [varName,var_struc.Type,var_struc.DimIds,natts] = netcdf.inqVar(ncid,varid-1);
    varName=matlab_sanitize_name (varName);
    if ~exist('v_names','var') || any(strcmp(varName,v_names))
      Vars.(varName) = var_struc;

      for attnum = 1:natts
        attName = netcdf.inqAttName(ncid,varid-1,attnum-1);
        [att_struc.Type,att_struc.Length] = netcdf.inqAtt(ncid,varid-1,attName);
        att_struc.Value = netcdf.getAtt(ncid,varid-1,attName);
        attName=matlab_sanitize_name (attName);

        Vars.(varName).Attr.(attName) = att_struc;
        if strcmp(attName,'FillValue'), Fill=att_struc.Value; end
        if strcmp(attName,'scale_factor'), scale=att_struc.Value; end
	if strcmp(attName,'add_offset'), offset=att_struc.Value; end
        if strcmp(attName,'missing_value'), Missing=att_struc.Value; end
      end
      datain = netcdf.getVar(ncid,varid-1);
      % At this point - should included a check for unsigned short or long integers
      %  In these cases - datain will wrap at 32767 and should have 65536 added
      %  The FillValue will be -1 and should be 65535
      % Need to know the xtype for unsigned integers
      if var_struc.Type == netcdf.getConstant('NC_CHAR')
          bad_data = [];
          missing_data = [];
      else
          bad_data = find(datain==Fill);
          missing_data = find(datain==Missing);
      end

      if scale~=1.; datain = double(datain) .* scale; end
      if offset~=0.; datain = double(datain) + double(offset); end
      if ~isempty(bad_data)
        datain=double(datain);
        datain(bad_data) = NaN;
      end
      if ~isempty(missing_data)
        datain=double(datain);
        datain(missing_data) = NaN;
      end
      Vars.(varName).Value = datain;
    end
  end
end

netcdf.close(ncid);

if ~exist('Vars','var'), Vars = []; end

function sanitized_name = matlab_sanitize_name (in_name)
%Create attribute and variable names that can be used as matlab field names

sanitized_name = in_name;

% Name can only contain letters, numbers or underscores
nonstr = find(~isletter(sanitized_name));
% Replace all invalid characters with underscores
if ~isempty(nonstr)
  for i=1:length(nonstr)
    if isnan(str2double(sanitized_name(nonstr(i))))
      sanitized_name(nonstr(i))='_';
    end
  end
end

% If the name has a leading underscore, remove it if possible
if (strcmp('_',sanitized_name(1)));
  if length(sanitized_name)>1, sanitized_name = sanitized_name(2:end); end
end

% If the first character is now not a letter, add a prefix.
% An attribute name of, say, '0h' is not permissible in matlab
if ~isletter(sanitized_name(1))
  sanitized_name = ['nc_' sanitized_name];
end

