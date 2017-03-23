function [f_memmap] = memmapfile_format % no input argument

n_hi = 20;  % No of hi rate data / block

f_memmap = { ...                             
    'int32',1,'utc_days'; ...      % UTC Time stamp of MDSR - days
    'uint32',1,'utc_secs'; ...     % UTC Time stamp of MDSR - seconds
    'uint32',1,'utc_usecs'; ...    % UTC Time stamp of MDSR - 10-6 secs
    'int16',1,'time_offset'; ...   % TAI - UTC offset
    'uint8',[2 1],'spare_char_1'; ... % Spare [2]
    'int32',[n_hi 1],'tdiff'; ...      % hi rate time difference of hi rate from 1Hz [n_hi]
    'int16',[n_hi 1],'tdiff_offset'; ...% hi rate TAI - UTC offset [n_hi]
    'uint32',1,'rec_cnt'; ...      % Record Counter
    'int32',1,'lat'; ...           % Geodetic Latitude (positive N)
    'int32',[n_hi 1],'hi_lat'; ...      % hi rate Geodetic Latitude (positive N) [n_hi]
    'int32',1,'lon'; ...            % Longitude (positive E, 0 at Greenwich)
    'int32',[n_hi 1],'hi_lon'; ...      % hi rate Longitude (positive E, 0 at Greenwich) [n_hi]
    'int32',1,'h'; ...             % Altitude of CoG above reference ellipsoid
    'int32',[n_hi 1],'hi_h'; ...       % hi rate altitude differences from 1 Hz altitude [n_hi].
    'int32',1,'h_rate'; ...        % Instantaneous altitude rate
    'uint32',[n_hi 1],'meas_conf_wd'; ... % Measurement Confidence Word
    'uint8',[2 1],'spare_char_2'; ...  % Spare [2]
    'int16',1,'pp'; ...            % Pulse Peakiness
    'int16',[n_hi 1],'hi_pp'; ...      % hi rate Pulse Peakiness [n_hi]
    'int16',[n_hi 1],'trk_error'; ...  % hi rate tracker error [n_hi]
    'uint32',1,'flag_trk'; ...     % hi rate tracker flag
    'uint8',[4 1],'spare_char_3'; ...  % Spare [4]
    ...% Range Measurements Group
    'uint32',1,'range'; ...             % ocean range
    'uint32',[n_hi 1],'hi_range'; ...       % hi rate ocean range [n_hi]
    'uint16',1,'range_sd'; ...          % Standard deviation of hi rate ocean range
    'uint16',1,'range_num'; ...         % Number of valid hi rate points for ocean range
    'uint32',1,'flag_range_av'; ...     % ocean range averaging flag
    'uint32',1,'range_ice'; ...         % ice ranges
    'uint32',[n_hi 1],'hi_range_ice'; ...   % hi rate ice2 ranges [n_hi]
    'uint16',1,'range_ice_sd'; ...      % Standard deviation of hi rate ice ranges
    'uint16',1,'range_ice_num'; ...     % Number of valid hi rate points for ice ranges
    'uint32',1,'flag_range_ice_av'; ... % ice range averaging flag
    ...% Range Corrections Group
    'int16',1,'corr_dopp'; ...          % Doppler correction
    'int16',1,'corr_uso'; ...           % range correction for USO drift
    'int16',1,'corr_cog'; ...        % range correction for Centre of Gravity
    'int16',1,'corr_int_cal1'; ...      % range correction for internal CAL1
    'int16',1,'corr_instr'; ...	        % range instrumental correction
    'int16',1,'corr_dry_trop'; ...	    % Model dry tropospheric correction
    'int16',1,'corr_wet_trop_mod'; ...  % Model wet tropospheric correction
    'int16',1,'corr_inv_baro'; ...      % Inverted barometer correction
    'int16',1,'corr_dac'; ...           % Dynamic Atmospheric Correction
    'int16',1,'corr_iono_gim'; ...      % Ionospheric correction from GIM model
    'int16',1,'corr_ssb'; ...           % Sea state bias
    'uint8',[6 1],'spare_char_4'; ...       % Spare [6]
    ...% SWH & Backscatter Measurements Group
    'int32',1,'swh_sq'; ...             % square of swh
    'int16',1,'swh'; ...                % Significant wave height
    'uint8',[2 1],'spare_char_5'; ...       % Spare [2]
    'int16',[n_hi 1],'hi_swh'; ...          % hi rate swh [n_hi]
    'uint16',1,'swh_sd'; ...            % Standard deviation of hi rate swh
    'uint16',1,'swh_num'; ...           % Number of valid hi rate points for swh
    'uint32',1,'flag_swh_av'; ...       % swh averaging flag
    'uint8',[2 1],'spare_char_6'; ...       % Spare [2]
    'int16',1,'sigma0'; ...                % corrected Ocean backscatter coefficient
    'int16',[n_hi 1],'hi_sigma0'; ...          % hi rate sigma0 [n_hi]
    'uint16',1,'sigma0_sd'; ...            % Standard deviation of hi rate sigma0
    'uint16',1,'sigma0_num'; ...           % Number of valid hi rate points for sigma0
    'uint32',1,'flag_sigma0_av'; ...       % sigma0 averaging flag
    'uint8',[2 1],'spare_char_7'; ...          % Spare [2]
    'int16',1,'sigma0_ice'; ...            % corrected Ice backscatter coefficient
    'int16',[n_hi 1],'hi_sigma0_ice'; ...      % hi rate ice sigma0 [n_hi]
    'uint16',1,'sigma0_ice_sd'; ...        % Standard deviation of hi rate ice sigma0
    'uint16',1,'sigma0_ice_num'; ...       % Number of valid hi rate points for ice sigma0
    'uint32',1,'flag_sigma0_ice_av'; ...   % ice sigma0 averaging flag
    'int32',1,'off_nad_ang_sq_wvform'; ... % Off nadir angle squared of the satellite from waveform data
    'uint8',[6 1],'spare_char_8'; ...          % Spare [6]
    'int16',1,'agc'; ...                   % AGC
    'int32',[n_hi 1],'hi_sigma0_scale'; ...    % hi rate scale factor for sigma0 [n_hi]
    'int16',1,'corr_swh_instr'; ...        % SWH instrument correction
    'int16',1,'corr_agc'; ...              % AGC correction
    'int16',1,'corr_sigma0_int_cal1'; ...  % internal cal1 correction on sigma0
    'int16',1,'corr_sigma0_instr'; ...     % instrument correction on sigma0
    'int16',1,'corr_atm_atten'; ...        % atmospheric attenuation correction on sigma0
    'uint8',[6 1],'spare_char_9'; ...          % Spare [6]
    ...% Geophysical Corrections Group
    'int32',1,'h_mss1'; ...                % Mean sea-surface height 1
    'int32',1,'h_mss2'; ...                % Mean sea-surface height 2
    'int32',1,'h_geoid'; ...               % Geoid height
    'int32',1,'h_dem'; ...                 % Ocean depth/land elevation (DEM)
    'int32',1,'h_mdt'; ...                 % Mean dynamic topography
    'uint8',[8 1],'spare_char_10'; ...         % Spare [8]
    'int16',1,'h_tot_geocen_ocn_tide_sol1'; ...	  % Total geocentric ocean tide height (solution 1)
    'int16',1,'h_tot_geocen_ocn_tide_sol2'; ...	  % Total geocentric ocean tide height (solution 2)
    'int16',1,'h_long_period_ocn_tide'; ...	      % Long period Tide height
    'int16',1,'h_long_period_ocn_tide_no_eq'; ... % Long period Tide height - none equlibrium tides
    'int16',1,'h_tide_load_sol1'; ...         % Tidal loading height (solution 1)
    'int16',1,'h_tide_load_sol2'; ...         % Tidal loading height (solution 2)
    'int16',1,'h_tide_solid'; ...             % Solid earth tide height
    'int16',1,'h_tide_geocen_pole'; ...       % Geocentric pole tide height
    'uint8',[6 1],'spare_char_11'; ...            % Spare [6]
    'int16',1,'wsp'; ...               % RA2 wind speed 
    'int16',1,'wsp_u_mod'; ...         % u component of the model wind vector
    'int16',1,'wsp_v_mod'; ...         % v component of the model wind vector 
    'uint16',1,'flag_surface_type'; ...% Surface type flag
    'uint8',[2 1],'spare_char_12' ...      % Spare [2]
};
return;
