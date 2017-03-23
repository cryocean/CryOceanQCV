is this used?
% read individual files for testing
clear all
data_type = 'SIR_GOP_L2';
dir_in=['/Volumes/scratch/general/cryosat/' data_type '/'];
% UPDATE HERE JAN 2017
% dir_in=['/noc/mpoc/cryosat/' data_type '/'];
dir_in=[dir_in,'2014/','04'];
fns = dir([dir_in,'/','CS*20140414*1.DBL']);
fns = struct2cell(fns);
fns = fns(1,:);
t2 = cell(length(fns),1);
for j=1:length(fns)
    fn = fns{j};
    
    mph_size = 1247;
    n_hi = 20;  % No of hi rate data / block
    infile = [dir_in '/' fn];
    
    f_id = fopen(infile,'r','ieee-be');
    if (f_id < 0)
        error('File Not Found');
    end
    
    % check if file is empty
    fseek(f_id, 0,'eof');
    if(ftell(f_id) == 0)
        DS = [];
        SPH = [];
        MPH = [];
        fclose(f_id);
        return
    end
    
    fseek(f_id, 0,'bof');
    
    % Read MPH - to character array
    fmt = sprintf('%%%dc',mph_size);
    MPH = fscanf(f_id,fmt,1);
    
    if(length(MPH) ~= mph_size)
        error('Length of MPH incorrect')
    end
    
    % Save useful parts of the MPH to l2_MPH structure
    MPH = esa_ph(MPH);
    
    fseek(f_id, mph_size,'bof');
    
    % Read SPH to character array - define format to read using MPH info
    fmt = sprintf('%%%dc',MPH.sph_size);
    SPH = fscanf(f_id,fmt,1);
    
    % Save useful parts of the SPH to l2_SPH structure
    SPH = esa_ph(SPH);
    
    i = 1;
    
    ds = c2_op_l2_by_var_fmc(SPH.DS(i).num_dsr);
    f_memmap = memmapfile_format;
    flds = fieldnames(ds);
    n_flds = length(flds);
    
    mem_map = memmapfile(infile,'Offset',SPH.DS(i).ds_offset, ...
        'Format',f_memmap,'Repeat',SPH.DS(i).num_dsr); % memory map
    mem_map = mem_map.Data;
    utc_time = (double(swapbytes([mem_map(:).utc_days]))*86400.) + ...
        double(swapbytes([mem_map(:).utc_secs]))+ ...
        (double(swapbytes([mem_map(:).utc_usecs]))*1e-6);
    time = utc_time+double(swapbytes([mem_map(:).time_offset]));
    hi_utctime = utc_time(ones(n_hi,1),:)+...
        double(swapbytes([mem_map(:).tdiff])).*1e-6;
    hi_time = hi_utctime+double(swapbytes([mem_map(:).tdiff_offset]));
    
    fclose(f_id);
    
    t = unique(time);
    if(length(t) ~= length(time))
        disp([length(t),length(time)])
        disp(fn)
    end
    t2{j} = time;
end


