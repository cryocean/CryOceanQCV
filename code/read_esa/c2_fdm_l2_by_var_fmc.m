function [rec] = c2_fdm_l2_by_var_fmc(n_recs)

% Define a structure to hold the SIR_L2_FDM records from a
% CRYOSAT L2 File, stored by variable
% [rec] = c2_op_l2_by_var(n_recs)

n_hi = 20;  % No of hi rate data / block

rec = struct ( ...
    ...% Time and orbit group
    'time', zeros(n_recs,1), ...                             % full TAI
    'tai_days', cast(zeros(n_recs,1), 'int32'), ...          % TAI Time stamp of MDSR - days
    'tai_secs', cast(zeros(n_recs,1), 'uint32'), ...         % TAI Time stamp of MDSR - seconds
    'tai_usecs', cast(zeros(n_recs,1), 'uint32'), ...        % TAI Time stamp of MDSR - 10-6 secs
    'utc_time', zeros(n_recs,1), ...                         % full UTC time
    'tdiff', cast(zeros(n_recs,n_hi), 'int32'), ...          % hi rate time difference of hi rate from 1Hz [n_hi]
    'hi_time', zeros(n_recs,n_hi), ...                       % hi rate full TAI [n_hi]
    'hi_utctime', zeros(n_recs,n_hi), ...                    % hi rate full UTC time [n_hi]
    'lat', cast(zeros(n_recs,1), 'int32'), ...               % Geodetic Latitude (positive N)
    'hi_lat', cast(zeros(n_recs,n_hi), 'int32'), ...         % hi rate Geodetic Latitude (positive N) [n_hi]
    'lon', cast(zeros(n_recs,1), 'int32'), ...               % Longitude (positive E, 0 at Greenwich)
    'hi_lon', cast(zeros(n_recs,n_hi), 'int32'), ...         % hi rate Longitude (positive E, 0 at Greenwich) [n_hi]
    'rec_cnt', cast(zeros(n_recs,1), 'uint32'), ...          % Record Counter
    'meas_conf_wd', cast(zeros(n_recs,1),'uint32'), ...      % Measurement Confidence Word
    'h', cast(zeros(n_recs,1),'int32'), ...                  % Altitude of CoG above reference ellipsoid
    'hi_h', cast(zeros(n_recs,n_hi),'int32'), ...            % 20Hz altitude of CoG [n_hi].
    'h_rate', cast(zeros(n_recs,1),'int16'), ...             % Altitude rate derived from orbit
    ...% Range Measurements Group
    'range', cast(zeros(n_recs,1),'uint32'), ...             % ocean range
    'hi_range', cast(zeros(n_recs,n_hi),'uint32'), ...       % hi rate ocean range [n_hi]
    'range_sd', cast(zeros(n_recs,1),'uint16'), ...          % Standard deviation of hi rate ocean range
    'range_num', cast(zeros(n_recs,1),'uint16'), ...         % Number of valid hi rate points for ocean range
    'flag_range_av', cast(zeros(n_recs,1),'uint32'), ...     % ocean range averaging flag
    'range_ocog', cast(zeros(n_recs,1),'uint32'), ...         % ocog ranges
    'hi_range_ocog', cast(zeros(n_recs,n_hi),'uint32'), ...   % hi rate ocog ranges [n_hi]
    'range_ocog_sd', cast(zeros(n_recs,1),'uint16'), ...      % Standard deviation of hi rate ocog ranges
    'range_ocog_num', cast(zeros(n_recs,1),'uint16'), ...     % Number of valid hi rate points for ocog ranges
    'flag_range_ocog_av', cast(zeros(n_recs,1),'uint32'), ... % ocog range averaging flag
    ...% Range Corrections Group
    'corr_dopp', cast(zeros(n_recs,1),'int16'), ...          % Doppler correction
    'corr_dry_trop', cast(zeros(n_recs,1),'int16'), ...	     % Model dry tropospheric correction
    'corr_wet_trop_mod', cast(zeros(n_recs,1),'int16'), ...  % Model wet tropospheric correction
    'corr_inv_baro', cast(zeros(n_recs,1),'int16'), ...      % Inverted barometer correction
    'corr_dac', cast(zeros(n_recs,1),'int16'), ...           % Dynamic Atmospheric Correction
    'corr_iono_gim', cast(zeros(n_recs,1),'int16'), ...      % Ionospheric correction from GIM model
    'corr_ssb', cast(zeros(n_recs,1),'int16'), ...           % Sea state bias
    ...% SWH & Backscatter Measurements Group
    'swh_sq', cast(zeros(n_recs,1),'int32'), ...             % square of swh
    'swh', cast(zeros(n_recs,1),'int16'), ...                % Significant wave height
    'hi_swh_sq', cast(zeros(n_recs,n_hi),'int32'), ...          % hi rate of squared swh [n_hi]
    'swh_sq_sd', cast(zeros(n_recs,1),'uint16'), ...            % Standard deviation of hi rate squared swh
    'swh_sq_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for squared swh
    'flag_swh_sq_av', cast(zeros(n_recs,1),'uint32'), ...       % squared swh averaging flag
    'sigma0', cast(zeros(n_recs,1),'int16'), ...                % corrected Ocean backscatter coefficient
    'hi_sigma0', cast(zeros(n_recs,n_hi),'int16'), ...          % hi rate sigma0 [n_hi]
    'sigma0_sd', cast(zeros(n_recs,1),'uint16'), ...            % Standard deviation of hi rate sigma0
    'sigma0_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for sigma0
    'flag_sigma0_av', cast(zeros(n_recs,1),'uint32'), ...       % sigma0 averaging flag
    'sigma0_ocog', cast(zeros(n_recs,1),'int16'), ...            % corrected ocog backscatter coefficient
    'hi_sigma0_ocog', cast(zeros(n_recs,n_hi),'int16'), ...      % hi rate ocog sigma0 [n_hi]
    'sigma0_ocog_sd', cast(zeros(n_recs,1),'uint16'), ...        % Standard deviation of hi rate ocog sigma0
    'sigma0_ocog_num', cast(zeros(n_recs,1),'uint16'), ...       % Number of valid hi rate points for ocog sigma0
    'flag_sigma0_ocog_av', cast(zeros(n_recs,1),'uint32'), ...   % ocog sigma0 averaging flag
    'off_nad_ang_platform', cast(zeros(n_recs,1),'int32'), ... % Off nadir angle of the satellite from platform data
    ...% Geophysical Corrections Group
    'h_mss1', cast(zeros(n_recs,1),'int32'), ...                % Mean sea-surface height 1
    'h_geoid', cast(zeros(n_recs,1),'int32'), ...               % Geoid height
    'h_dem', cast(zeros(n_recs,1),'int32'), ...                 % Ocean depth/land elevation (DEM)
    'h_tot_geocen_ocn_tide_sol1', cast(zeros(n_recs,1),'int16'), ...	% Total geocentric ocean tide height
    'h_long_period_ocn_tide', cast(zeros(n_recs,1),'int16'), ...	    % Long period Tide height
    'h_tide_load_sol1', cast(zeros(n_recs,1),'int16'), ...         % Tidal loading height
    'h_tide_solid', cast(zeros(n_recs,1),'int16'), ...             % Solid earth tide height
    'h_tide_geocen_pole', cast(zeros(n_recs,1),'int16'), ...       % Geocentric pole tide height
    'wsp', cast(zeros(n_recs,1),'int16'), ...               % RA2 wind speed 
    'wsp_u_mod', cast(zeros(n_recs,1),'int16'), ...         % u component of the model wind vector
    'wsp_v_mod', cast(zeros(n_recs,1),'int16'), ...         % v component of the model wind vector 
    'hi_pp', cast(zeros(n_recs,n_hi),'int16'), ...           % hi rate Ku-band Peakiness [n_hi]
    'retrk_flag', cast(zeros(n_recs,1),'uint32'), ...         % ocean retracking quality flag
    'flag_surface_type', cast(zeros(n_recs,1),'uint16') ... % Surface type flag
);
