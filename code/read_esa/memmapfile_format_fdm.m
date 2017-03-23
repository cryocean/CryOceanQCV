function [f_memmap] = memmapfile_format_fdm % no input argument

n_hi = 20;  % No of hi rate data / block

f_memmap = { ...                             
    'int32',1,'tai_days'; ...      % TAI Time stamp of MDSR - days
    'uint32',1,'tai_secs'; ...     % TAI Time stamp of MDSR - seconds
    'uint32',1,'tai_usecs'; ...    % TAI Time stamp of MDSR - 10-6 secs
    'int32',[n_hi 1],'tdiff'; ...      % hi rate time difference of hi rate from 1Hz [n_hi]
    'int32',1,'lat'; ...           % Geodetic Latitude (positive N)
    'int32',[n_hi 1],'hi_lat'; ...      % hi rate Geodetic Latitude (positive N) [n_hi]
    'int32',1,'lon'; ...            % Longitude (positive E, 0 at Greenwich)
    'int32',[n_hi 1],'hi_lon'; ...      % hi rate Longitude (positive E, 0 at Greenwich) [n_hi]
    'uint32',1,'rec_cnt'; ...      % Record Counter
    'uint32',1,'meas_conf_wd'; ... % Measurement Confidence Word
    'int32',1,'h'; ...             % Altitude of CoG above reference ellipsoid
    'int32',[n_hi 1],'hi_h'; ...       % hi rate altitude differences from 1 Hz altitude [n_hi].
    'int16',1,'h_rate'; ...        % Instantaneous altitude rate
    'uint16',1,'spare_char_2'; ...  % Spare
    ...% Range Measurements Group
    'uint32',1,'range'; ...             % ocean range
    'uint32',[n_hi 1],'hi_range'; ...       % hi rate ocean range [n_hi]
    'uint16',1,'range_sd'; ...          % Standard deviation of hi rate ocean range
    'uint16',1,'range_num'; ...         % Number of valid hi rate points for ocean range
    'uint32',1,'flag_range_av'; ...     % ocean range averaging flag
    'uint32',1,'range_ocog'; ...         % ocog ranges
    'uint32',[n_hi 1],'hi_range_ocog'; ...   % hi rate ocog ranges [n_hi]
    'uint16',1,'range_ocog_sd'; ...      % Standard deviation of hi rate ocog ranges
    'uint16',1,'range_ocog_num'; ...     % Number of valid hi rate points for ocog ranges
    'uint32',1,'flag_range_ocog_av'; ... % ocog range averaging flag
    ...% Range Corrections Group
    'int16',1,'corr_dopp'; ...          % Doppler correction
    'int16',1,'corr_dry_trop'; ...	    % Model dry tropospheric correction
    'int16',1,'corr_wet_trop_mod'; ...  % Model wet tropospheric correction
    'int16',1,'corr_inv_baro'; ...      % Inverted barometer correction
    'int16',1,'corr_dac'; ...           % Dynamic Atmospheric Correction
    'int16',1,'corr_iono_gim'; ...      % Ionospheric correction from GIM model
    'int16',1,'corr_ssb'; ...           % Sea state bias
    'uint16',[3 1],'spare_char_4'; ...       % Spare [6]
    ...% SWH & Backscatter Measurements Group
    'int32',1,'swh_sq'; ...             % square of swh
    'int16',1,'swh'; ...                % Significant wave height
    'uint16',1,'spare_char_5'; ...       % Spare
    'int32',[n_hi 1],'hi_swh_sq'; ...          % hi rate swh [n_hi]
    'uint16',1,'swh_sq_sd'; ...            % Standard deviation of hi rate squared swh
    'uint16',1,'swh_sq_num'; ...           % Number of valid hi rate points for squared swh
    'uint32',1,'flag_swh_sq_av'; ...       % Squared swh averaging flag
    'uint16',1,'spare_char_6'; ...       % Spare
    'int16',1,'sigma0'; ...                % corrected Ocean backscatter coefficient
    'int16',[n_hi 1],'hi_sigma0'; ...          % hi rate sigma0 [n_hi]
    'uint16',1,'sigma0_sd'; ...            % Standard deviation of hi rate sigma0
    'uint16',1,'sigma0_num'; ...           % Number of valid hi rate points for sigma0
    'uint32',1,'flag_sigma0_av'; ...       % sigma0 averaging flag
    'uint16',1,'spare_char_7'; ...          % Spare
    'int16',1,'sigma0_ocog'; ...            % corrected ocog backscatter coefficient
    'int16',[n_hi 1],'hi_sigma0_ocog'; ...      % hi rate ocog sigma0 [n_hi]
    'uint16',1,'sigma0_ocog_sd'; ...        % Standard deviation of hi rate ocog sigma0
    'uint16',1,'sigma0_ocog_num'; ...       % Number of valid hi rate points for ocog sigma0
    'uint32',1,'flag_sigma0_ocog_av'; ...   % ocog sigma0 averaging flag
    'int32',1,'off_nad_ang_platform'; ... % Off nadir angle of the satellite from platform data
    'uint32',1,'spare_char_8'; ...          % Spare
    ...% Geophysical Corrections Group
    'int32',1,'h_mss1'; ...                % Mean sea-surface height 1
    'int32',1,'h_geoid'; ...               % Geoid height
    'int32',1,'h_dem'; ...                 % Ocean depth/land elevation (DEM)
    'int16',1,'h_tot_geocen_ocn_tide_sol1'; ...	  % Total geocentric ocean tide height
    'int16',1,'h_long_period_ocn_tide'; ...	      % Long period Tide height
    'int16',1,'h_tide_load_sol1'; ...         % Tidal loading height
    'int16',1,'h_tide_solid'; ...             % Solid earth tide height
    'int16',1,'h_tide_geocen_pole'; ...       % Geocentric pole tide height
    'int16',1,'wsp'; ...               % RA2 wind speed 
    'int16',1,'wsp_u_mod'; ...         % u component of the model wind vector
    'int16',1,'wsp_v_mod'; ...         % v component of the model wind vector 
    'int16',[n_hi,1],'hi_pp'; ...
    'uint32',1,'retrk_flag'; ...
    'uint16',1,'flag_surface_type'; ...% Surface type flag
    'uint16',1,'spare_char_12' ...      % Spare
};
return;
