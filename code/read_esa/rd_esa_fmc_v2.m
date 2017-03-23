function [DS,SPH,MPH] = rd_esa_fmc_v2(dir_in,fn,scaled)
%%Read ALL measurement datasets from ESA Cryosat files
%
% USAGE : [DS,SPH,MPH] = rd_esa_fmc_v2(dir_in,fn,scaled)
%
%
% INPUT :    dir_in : path to binary data
%            fn : file name
%            scaled :  integer of value 1 (scale) or 0 (do not scale)
%
% Note1: This function assumes that the byte ordering of the data to read is
% different from that of your operating system (hence the use of "swapbytes").
% If the data and your OS have the same endinness then the use of swapbytes
% is not necessary.
%
% Note2: For FDM, the conversion from TAU to UTC is done assuming a
% difference of 35s. At present, if a new leap second is inserted this
% needs to be modified manually since there is no information in the FDM
% files

%% Error check input arguments and set verbosity
if (nargin < 2)
    error('Too few input arguments');
end

%% Open file
mph_size = 1247;
n_hi = 20;  % No of hi rate data / block
infile = [dir_in '/' fn];
fprintf(1,'File: %s\n',infile);

f_id = fopen(infile,'r','ieee-be');
if (f_id < 0)
    error('File Not Found');
end

%% check if file is empty
fseek(f_id, 0,'eof'); % go to end of file
if(ftell(f_id) == 0) % position in file is 1 (ie file empty)/
   DS = [];
   SPH = [];
   MPH = [];
   fclose(f_id);
   return
end

%% Read MPH
% Move to start of MPH
fseek(f_id, 0,'bof');

% Read MPH - to character array
fmt = sprintf('%%%dc',mph_size);
MPH = fscanf(f_id,fmt,1);

if(length(MPH) ~= mph_size)
    error('Length of MPH incorrect')
end

% Save useful parts of the MPH to l2_MPH structure
MPH = esa_ph(MPH);

%% Read SPH
% Move to start of SPH
fseek(f_id, mph_size,'bof');

% Read SPH to character array - define format to read using MPH info
fmt = sprintf('%%%dc',MPH.sph_size);
SPH = fscanf(f_id,fmt,1);

% Save useful parts of the SPH to l2_SPH structure
SPH = esa_ph(SPH);


%% Read all the datasets
for i=1:length(SPH.DS)
    %% For Measurement DataSets:- read and save variables
    if strcmp(SPH.DS(i).ds_type,'M')
        
        % Pre-define Data structure: formatted by variable not time
        switch SPH.DS(i).ds_name
            case {'SIR_L2_IOP','SIR_L2_GOP'}
                ds = c2_op_l2_by_var_fmc(SPH.DS(i).num_dsr);
                f_memmap = memmapfile_format;
            case {'SIR_FDM_L2'}
                ds = c2_fdm_l2_by_var_fmc(SPH.DS(i).num_dsr);
                f_memmap = memmapfile_format_fdm;
            otherwise
                display(['Sorry, cannot read this DataType ' SPH.DS(i).ds_name]);
                return
        end
        
        % Find fields to work with for this MDSR
        flds = fieldnames(ds);
        n_flds = length(flds);
        
        mem_map = memmapfile(infile,'Offset',SPH.DS(i).ds_offset, ...
            'Format',f_memmap,'Repeat',SPH.DS(i).num_dsr); % memory map
        mem_map = mem_map.Data;
        
        switch SPH.DS(i).ds_name
            case {'SIR_L2_IOP','SIR_L2_GOP'}
                % compute time
                utc_time = (double(swapbytes([mem_map(:).utc_days]))*86400.) + ...
                    double(swapbytes([mem_map(:).utc_secs]))+ ...
                    (double(swapbytes([mem_map(:).utc_usecs]))*1e-6);
                time = utc_time+double(swapbytes([mem_map(:).time_offset]));
                hi_utctime = utc_time(ones(n_hi,1),:)+...
                    double(swapbytes([mem_map(:).tdiff])).*1e-6;
                hi_time = hi_utctime+double(swapbytes([mem_map(:).tdiff_offset]));
            case {'SIR_FDM_L2'}
                % compute time
                time = (double(swapbytes([mem_map(:).tai_days]))*86400.) + ...
                    double(swapbytes([mem_map(:).tai_secs]))+ ...
                    (double(swapbytes([mem_map(:).tai_usecs])).*1e-6);
                utc_time = time-36; % TAI is currently ahead of UTC by 36 seconds
                hi_time = time(ones(n_hi,1),:)+...
                    double(swapbytes([mem_map(:).tdiff])).*1e-6;
                hi_utctime = hi_time-36;
        end
        time = num2cell(time,1);
        utc_time = num2cell(utc_time,1);
        hi_time = num2cell(hi_time,1);
        hi_utctime = num2cell(hi_utctime,1);
        
        % store the computed time
        [mem_map(:).time] = time{:}; % full TAI time
        [mem_map(:).utc_time] = utc_time{:}; % full UTC time
        [mem_map(:).hi_time] = hi_time{:};
        [mem_map(:).hi_utctime] = hi_utctime{:};
        
        % Save useful parts of the record to the 'by variable' structure
        for k=1:n_flds
            fld = char(flds(k));
            if(any(strcmp(fld,{'time','utc_time','hi_time','hi_utctime'})))
                if (size(ds.(fld),1)>1)
                    ds.(fld)(:,:) = [mem_map(:).(fld)]';
                else
                    ds.(fld)(:) = [mem_map(:).(fld)];
                end
            else
                if (size(ds.(fld),1)>1)
                    ds.(fld)(:,:) = swapbytes([mem_map(:).(fld)]');
                else
                    ds.(fld)(:) = swapbytes([mem_map(:).(fld)]);
                end
            end
        end
        
        % Generate scaled variables if necessary
        if scaled
            switch SPH.DS(i).ds_name
                case {'SIR_L2_IOP','SIR_L2_GOP'}
                    ds = c2_op_l2_scaled(ds,SPH.DS(i).num_dsr);
                case {'SIR_FDM_L2'}
                    ds = c2_fdm_l2_scaled(ds,SPH.DS(i).num_dsr);
            end
        end

        %Save to DS structure
        DS.(SPH.DS(i).ds_name) = ds;
        
    end
    
end

%% Close file
fclose(f_id);
