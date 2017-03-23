function [rec] = c2_op_l2_by_var_fmc(n_recs)

% Define a structure to hold the SIR_L2_IOP / SIR_L2_GOP records from a
% CRYOSAT L2 (I/G)OP File, stored by variable
% [rec] = c2_op_l2_by_var(n_recs)
%
% For high rate (n_hi) bit values:
% First n_hi least significant bits (bits 0-n_hi-1) correspond to the n_hi values
% (one per data block) containing:
% 0=valid, 1=invalid.
% Unused bits set to 1.
% Bit 0 applies to the first data block.
n_hi = 20;  % No of hi rate data / block

rec = struct ( ...
    ...% Time and orbit group
    'time', zeros(n_recs,1), ...                             % full TAI
    'utc_days', cast(zeros(n_recs,1), 'int32'), ...          % UTC Time stamp of MDSR - days
    'utc_secs', cast(zeros(n_recs,1), 'uint32'), ...         % UTC Time stamp of MDSR - seconds
    'utc_usecs', cast(zeros(n_recs,1), 'uint32'), ...        % UTC Time stamp of MDSR - 10-6 secs
    'time_offset', cast(zeros(n_recs,1), 'int16'), ...       % TAI - UTC offset
    'utc_time', zeros(n_recs,1), ...                         % full UTC time
    'tdiff', cast(zeros(n_recs,n_hi), 'int32'), ...          % hi rate time difference of hi rate from 1Hz [n_hi]
    'hi_time', zeros(n_recs,n_hi), ...                       % hi rate full TAI [n_hi]
    'tdiff_offset', cast(zeros(n_recs,n_hi), 'int16'), ...   % hi rate TAI - UTC offset [n_hi]
    'hi_utctime', zeros(n_recs,n_hi), ...                    % hi rate full UTC time [n_hi]
    'rec_cnt', cast(zeros(n_recs,1), 'uint32'), ...          % Record Counter
    'lat', cast(zeros(n_recs,1), 'int32'), ...               % Geodetic Latitude (positive N)
    'hi_lat', cast(zeros(n_recs,n_hi), 'int32'), ...         % hi rate Geodetic Latitude (positive N) [n_hi]
    'lon', cast(zeros(n_recs,1), 'int32'), ...               % Longitude (positive E, 0 at Greenwich)
    'hi_lon', cast(zeros(n_recs,n_hi), 'int32'), ...         % hi rate Longitude (positive E, 0 at Greenwich) [n_hi]
    'h', cast(zeros(n_recs,1),'int32'), ...                  % Altitude of CoG above reference ellipsoid
    'hi_h', cast(zeros(n_recs,n_hi),'int32'), ...            % hi rate altitude differences from 1 Hz altitude [n_hi].
    'h_rate', cast(zeros(n_recs,1),'int32'), ...             % Instantaneous altitude rate
    'meas_conf_wd', cast(zeros(n_recs,n_hi),'uint32'), ...   % Measurement Confidence Word
    'pp', cast(zeros(n_recs,1),'int16'), ...                 % Pulse Peakiness
    'hi_pp', cast(zeros(n_recs,n_hi),'int16'), ...           % hi rate Pulse Peakiness [n_hi]
    'trk_error', cast(zeros(n_recs,n_hi),'int16'), ...       % hi rate tracker error [n_hi]
    'flag_trk', cast(zeros(n_recs,1),'uint32'), ...          % hi rate tracker flag [n_hi]
    ...% Range Measurements Group
    'range', cast(zeros(n_recs,1),'uint32'), ...             % ocean range
    'hi_range', cast(zeros(n_recs,n_hi),'uint32'), ...       % hi rate ocean range [n_hi]
    'range_sd', cast(zeros(n_recs,1),'uint16'), ...          % Standard deviation of hi rate ocean range
    'range_num', cast(zeros(n_recs,1),'uint16'), ...         % Number of valid hi rate points for ocean range
    'flag_range_av', cast(zeros(n_recs,1),'uint32'), ...     % ocean range averaging flag
    'range_ice', cast(zeros(n_recs,1),'uint32'), ...         % ice ranges
    'hi_range_ice', cast(zeros(n_recs,n_hi),'uint32'), ...   % hi rate ice2 ranges [n_hi]
    'range_ice_sd', cast(zeros(n_recs,1),'uint16'), ...      % Standard deviation of hi rate ice ranges
    'range_ice_num', cast(zeros(n_recs,1),'uint16'), ...     % Number of valid hi rate points for ice ranges
    'flag_range_ice_av', cast(zeros(n_recs,1),'uint32'), ... % ice range averaging flag
    ...% Range Corrections Group
    'corr_dopp', cast(zeros(n_recs,1),'int16'), ...          % Doppler correction
    'corr_uso', cast(zeros(n_recs,1),'int16'), ...           % range correction for USO drift
    'corr_cog', cast(zeros(n_recs,1),'int16'), ...           % range correction for Centre of Gravity
    'corr_int_cal1', cast(zeros(n_recs,1),'int16'), ...      % range correction for internal CAL1
    'corr_instr', cast(zeros(n_recs,1),'int16'), ...	     % range instrumental correction
    'corr_dry_trop', cast(zeros(n_recs,1),'int16'), ...	     % Model dry tropospheric correction
    'corr_wet_trop_mod', cast(zeros(n_recs,1),'int16'), ...  % Model wet tropospheric correction
    'corr_inv_baro', cast(zeros(n_recs,1),'int16'), ...      % Inverted barometer correction
    'corr_dac', cast(zeros(n_recs,1),'int16'), ...           % Dynamic Atmospheric Correction
    'corr_iono_gim', cast(zeros(n_recs,1),'int16'), ...      % Ionospheric correction from GIM model
    'corr_ssb', cast(zeros(n_recs,1),'int16'), ...           % Sea state bias
    ...% SWH & Backscatter Measurements Group
    'swh_sq', cast(zeros(n_recs,1),'int32'), ...             % square of swh
    'swh', cast(zeros(n_recs,1),'int16'), ...                % Significant wave height
    'hi_swh', cast(zeros(n_recs,n_hi),'int16'), ...          % hi rate swh [n_hi]
    'swh_sd', cast(zeros(n_recs,1),'uint16'), ...            % Standard deviation of hi rate swh
    'swh_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for swh
    'flag_swh_av', cast(zeros(n_recs,1),'uint32'), ...       % swh averaging flag
    'sigma0', cast(zeros(n_recs,1),'int16'), ...                % corrected Ocean backscatter coefficient
    'hi_sigma0', cast(zeros(n_recs,n_hi),'int16'), ...          % hi rate sigma0 [n_hi]
    'sigma0_sd', cast(zeros(n_recs,1),'uint16'), ...            % Standard deviation of hi rate sigma0
    'sigma0_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for sigma0
    'flag_sigma0_av', cast(zeros(n_recs,1),'uint32'), ...       % sigma0 averaging flag
    'sigma0_ice', cast(zeros(n_recs,1),'int16'), ...            % corrected Ice backscatter coefficient
    'hi_sigma0_ice', cast(zeros(n_recs,n_hi),'int16'), ...      % hi rate ice sigma0 [n_hi]
    'sigma0_ice_sd', cast(zeros(n_recs,1),'uint16'), ...        % Standard deviation of hi rate ice sigma0
    'sigma0_ice_num', cast(zeros(n_recs,1),'uint16'), ...       % Number of valid hi rate points for ice sigma0
    'flag_sigma0_ice_av', cast(zeros(n_recs,1),'uint32'), ...   % ice sigma0 averaging flag
    'off_nad_ang_sq_wvform', cast(zeros(n_recs,1),'int32'), ... % Off nadir angle squared of the satellite from waveform data
    'agc', cast(zeros(n_recs,1),'int16'), ...                   % AGC
    'hi_sigma0_scale', cast(zeros(n_recs,n_hi),'int32'), ...    % hi rate scale factor for sigma0 [n_hi]
    'corr_swh_instr', cast(zeros(n_recs,1),'int16'), ...        % SWH instrument correction
    'corr_agc', cast(zeros(n_recs,1),'int16'), ...              % AGC correction
    'corr_sigma0_int_cal1', cast(zeros(n_recs,1),'int16'), ...  % internal cal1 correction on sigma0
    'corr_sigma0_instr', cast(zeros(n_recs,1),'int16'), ...     % instrument correction on sigma0
    'corr_atm_atten', cast(zeros(n_recs,1),'int16'), ...        % atmospheric attenuation correction on sigma0
    ...% Geophysical Corrections Group
    'h_mss1', cast(zeros(n_recs,1),'int32'), ...                % Mean sea-surface height 1
    'h_mss2', cast(zeros(n_recs,1),'int32'), ...                % Mean sea-surface height 2
    'h_geoid', cast(zeros(n_recs,1),'int32'), ...               % Geoid height
    'h_dem', cast(zeros(n_recs,1),'int32'), ...                 % Ocean depth/land elevation (DEM)
    'h_mdt', cast(zeros(n_recs,1),'int32'), ...                 % Mean dynamic topography
    'h_tot_geocen_ocn_tide_sol1', cast(zeros(n_recs,1),'int16'), ...	% Total geocentric ocean tide height (solution 1)
    'h_tot_geocen_ocn_tide_sol2', cast(zeros(n_recs,1),'int16'), ...	% Total geocentric ocean tide height (solution 2)
    'h_long_period_ocn_tide', cast(zeros(n_recs,1),'int16'), ...	    % Long period Tide height
    'h_long_period_ocn_tide_no_eq', cast(zeros(n_recs,1),'int16'), ...	% Long period Tide height - none equlibrium tides
    'h_tide_load_sol1', cast(zeros(n_recs,1),'int16'), ...         % Tidal loading height (solution 1)
    'h_tide_load_sol2', cast(zeros(n_recs,1),'int16'), ...         % Tidal loading height (solution 2)
    'h_tide_solid', cast(zeros(n_recs,1),'int16'), ...             % Solid earth tide height
    'h_tide_geocen_pole', cast(zeros(n_recs,1),'int16'), ...       % Geocentric pole tide height
    'wsp', cast(zeros(n_recs,1),'int16'), ...               % RA2 wind speed 
    'wsp_u_mod', cast(zeros(n_recs,1),'int16'), ...         % u component of the model wind vector
    'wsp_v_mod', cast(zeros(n_recs,1),'int16'), ...         % v component of the model wind vector 
    'flag_surface_type', cast(zeros(n_recs,1),'uint16') ... % Surface type flag
);
