function [v] = c2_fdm_l2_scaled(v_in,n_recs)

% Scale the SIR_L2_FDM records from a
%  CRYOSAT L2 (I/G)OP File, stored by variable,
%  to SI units, setting default values to NaNs
% 'map' flags converted to arrays
% [v] = c2_fdm_l2_scaled(v_in,n_recs)
%
% For high rate (n_hi) bit values:
% First n_hi least significant bits (bits 0-n_hi-1) correspond to the n_hi values
% (one per data block) containing:
% 0=valid, 1=invalid.
% Unused bits set to 1.
% Bit 0 applies to the first data block.

%% setup structure to hold scaled variables
n_hi = 20;  % No of hi rate data / block

v = struct ( ...
    ...% Time and orbit group
    'time', zeros(n_recs,1), ...                             % full TAI
    'utc_time', zeros(n_recs,1), ...                         % full UTC time
    'tdiff', zeros(n_recs,n_hi), ...          % hi rate time difference of hi rate from 1Hz [n_hi]
    'hi_time', zeros(n_recs,n_hi), ...                       % hi rate full TAI [n_hi]
    'hi_utctime', zeros(n_recs,n_hi), ...                    % hi rate full UTC time [n_hi]
    'lat', zeros(n_recs,1), ...               % Geodetic Latitude (positive N)
    'hi_lat',zeros(n_recs,n_hi), ...         % hi rate Geodetic Latitude (positive N) [n_hi]
    'lon', zeros(n_recs,1), ...               % Longitude (positive E, 0 at Greenwich)
    'hi_lon', zeros(n_recs,n_hi), ...         % hi rate Longitude (positive E, 0 at Greenwich) [n_hi]
    'rec_cnt', cast(zeros(n_recs,1), 'uint32'), ...          % Record Counter
    'meas_conf_wd', cast(zeros(n_recs,1),'uint32'), ...      % Measurement Confidence Word
    'h', zeros(n_recs,1), ...                  % Altitude of CoG above reference ellipsoid
    'hi_h', zeros(n_recs,n_hi), ...            % 20Hz altittude of CoG [n_hi].
    'h_rate', zeros(n_recs,1), ...             % Instantaneous altitude rate
    'flag_mcd', struct( ...
          'block_degraded', cast(zeros(n_recs,1),'int8'), ...
          'blank_block', cast(zeros(n_recs,1),'int8'), ...
          'datation_degrad',cast(zeros(n_recs,1),'int8'), ...
          'orbit_propag_err', cast(zeros(n_recs,1),'int8'), ...
          'orbit_file_change', cast(zeros(n_recs,1),'int8'), ...
          'orbit_discontinuity', cast(zeros(n_recs,1),'int8'), ...
          'echo_saturation', cast(zeros(n_recs,1),'int8'), ...
          'other_echo_err', cast(zeros(n_recs,1),'int8'), ...
          'rx_ch1_err',cast(zeros(n_recs,1),'int8'), ...
          'rx_ch2_err',cast(zeros(n_recs,1),'int8'), ...
          'window_delay_incon',cast(zeros(n_recs,1),'int8'), ...
          'agc_incon',cast(zeros(n_recs,1),'int8'), ...
          'cal1_corr_miss', cast(zeros(n_recs,1),'int8'), ...
          'cal1_corr_ipf', cast(zeros(n_recs,1),'int8'), ...
          'doris_uso_corr', cast(zeros(n_recs,1),'int8'), ...
          'complex_call_ipf',cast(zeros(n_recs,1),'int8'), ...
          'trk_echo_err', cast(zeros(n_recs,1),'int8'), ...
          'rx1_echo_err', cast(zeros(n_recs,1),'int8'), ...
          'rx2_echo_err', cast(zeros(n_recs,1),'int8'), ...
          'noise_pow_incon', cast(zeros(n_recs,1),'int8'), ...
          'azim_cal_miss', cast(zeros(n_recs,1),'int8'), ...
          'default_azim_cal_ipf', cast(zeros(n_recs,1),'int8'), ...
          'range_cal_miss', cast(zeros(n_recs,1),'int8'), ...
          'range_cal_ipf', cast(zeros(n_recs,1),'int8'), ...
          'cal2_corr_miss', cast(zeros(n_recs,1),'int8'), ...
          'cal2_corr_ipf', cast(zeros(n_recs,1),'int8'), ...
          'power_scaling_err', cast(zeros(n_recs,1),'int8'), ...
          'attitude_corr_miss', cast(zeros(n_recs,1),'int8'), ...
          'attitude_corr_err', cast(zeros(n_recs,1),'int8') ...
         ), ...
    ...% Range Measurements Group
    'range', zeros(n_recs,1), ...             % ocean range
    'hi_range', zeros(n_recs,n_hi), ...       % hi rate ocean range [n_hi]
    'range_sd', zeros(n_recs,1), ...          % Standard deviation of hi rate ocean range
    'range_num', cast(zeros(n_recs,1),'uint16'), ...         % Number of valid hi rate points for ocean range
    'flag_range_av', cast(zeros(n_recs,n_hi+1),'uint8'), ...     % ocean range averaging flag
    'range_ocog', zeros(n_recs,1), ...         % ocog ranges
    'hi_range_ocog', zeros(n_recs,n_hi), ...   % hi rate ocog ranges [n_hi]
    'range_ocog_sd', zeros(n_recs,1), ...      % Standard deviation of hi rate ocog ranges
    'range_ocog_num', cast(zeros(n_recs,1),'uint16'), ...     % Number of valid hi rate points for ocog ranges
    'flag_range_ocog_av', cast(zeros(n_recs,n_hi+1),'uint8'), ... % ocog range averaging flag
    ...% Range Corrections Group
    'corr_dopp', cast(zeros(n_recs,1),'int16'), ...          % Doppler correction
    'corr_dry_trop', cast(zeros(n_recs,1),'int16'), ...	     % Model dry tropospheric correction
    'corr_wet_trop_mod', cast(zeros(n_recs,1),'int16'), ...  % Model wet tropospheric correction
    'corr_inv_baro', cast(zeros(n_recs,1),'int16'), ...      % Inverted barometer correction
    'corr_dac', cast(zeros(n_recs,1),'int16'), ...           % Dynamic Atmospheric Correction
    'corr_iono_gim', cast(zeros(n_recs,1),'int16'), ...      % Ionospheric correction from GIM model
    'corr_ssb', cast(zeros(n_recs,1),'int16'), ...           % Sea state bias
    ...% SWH & Backscatter Measurements Group
    'swh_sq',zeros(n_recs,1), ...             % square of swh
    'swh', zeros(n_recs,1), ...                % Significant wave height
    'hi_swh_sq', zeros(n_recs,n_hi), ...          % hi rate squared swh [n_hi]
    'swh_sq_sd', zeros(n_recs,1), ...            % Standard deviation of hi rate squared swh
    'swh_sq_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for squared swh
    'flag_swh_sq_av', cast(zeros(n_recs,n_hi+1),'uint8'), ...       % Squared swh averaging flag
    'sigma0', zeros(n_recs,1), ...                % corrected Ocean backscatter coefficient
    'hi_sigma0', zeros(n_recs,n_hi), ...          % hi rate sigma0 [n_hi]
    'sigma0_sd', zeros(n_recs,1), ...            % Standard deviation of hi rate sigma0
    'sigma0_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for sigma0
    'flag_sigma0_av', cast(zeros(n_recs,n_hi+1),'uint8'), ...       % sigma0 averaging flag
    'sigma0_ocog', zeros(n_recs,1), ...            % corrected ocog backscatter coefficient
    'hi_sigma0_ocog', zeros(n_recs,n_hi), ...      % hi rate ocog sigma0 [n_hi]
    'sigma0_ocog_sd', zeros(n_recs,1), ...        % Standard deviation of hi rate ocog sigma0
    'sigma0_ocog_num', cast(zeros(n_recs,1),'uint16'), ...       % Number of valid hi rate points for ocog sigma0
    'flag_sigma0_ocog_av', cast(zeros(n_recs,n_hi+1),'uint8'), ...   % ocog sigma0 averaging flag
    'off_nad_ang_platform', zeros(n_recs,1), ... % Off nadir angle of the satellite from platform data
    ...% Geophysical Corrections Group
    'h_mss1',zeros(n_recs,1), ...                % Mean sea-surface height 1
    'h_geoid', zeros(n_recs,1), ...               % Geoid height
    'h_dem', zeros(n_recs,1), ...                 % Ocean depth/land elevation (DEM)
    'h_tot_geocen_ocn_tide_sol1', zeros(n_recs,1), ...	% Total geocentric ocean tide height
    'h_long_period_ocn_tide', zeros(n_recs,1), ...	    % Long period Tide height
    'h_tide_load_sol1', zeros(n_recs,1), ...         % Tidal loading height
    'h_tide_solid', zeros(n_recs,1), ...             % Solid earth tide height
    'h_tide_geocen_pole', zeros(n_recs,1), ...       % Geocentric pole tide height
    'wsp', zeros(n_recs,1), ...               % RA2 wind speed 
    'wsp_u_mod', zeros(n_recs,1), ...         % u component of the model wind vector
    'wsp_v_mod', zeros(n_recs,1), ...         % v component of the model wind vector 
    'hi_pp', zeros(n_recs,n_hi), ...           % hi rate Ku-band Peakiness [n_hi]
    'retrk_flag', cast(zeros(n_recs,n_hi),'uint8'), ...         % ocean retracking quality flag
    'flag_surface_type', cast(zeros(n_recs,1),'uint16') ... % Surface type flag
);

%% Define which variables are left as they are, scaled or mapped
noscale={'time','tdiff','utc_time',...
          'hi_time','hi_utctime', ...
          'rec_cnt','meas_conf_wd', ...
          'range_num','range_ocog_num', ...
          'swh_sq_num','sigma0_num', ...
          'sigma0_ocog_num','hi_pp', 'flag_surface_type',...
         };
notflag={'rec_cnt','range_num','range_ocog_num' ...
         'swh_sq_num','sigma0_num','sigma0_ocog_num' ...
         };
scale01={};
scale02={
         'sigma0','hi_sigma0','sigma0_sd', ...
         'sigma0_ocog','hi_sigma0_ocog','sigma0_ocog_sd'...
         };
scale03={'h','hi_h','h_rate',...
         'range','hi_range','range_sd', ...
         'range_ocog','hi_range_ocog','range_ocog_sd', ...
         'corr_dopp','corr_dry_trop','corr_wet_trop_mod', ...
         'corr_inv_baro','corr_dac', ...
         'corr_iono_gim','corr_ssb', ...
         'swh','h_mss1','h_geoid','h_dem', ...
         'h_tot_geocen_ocn_tide_sol1', ...
         'h_long_period_ocn_tide', ...
         'h_tide_load_sol1',...
         'h_tide_solid','h_tide_geocen_pole',...
         'wsp','wsp_u_mod','wsp_v_mod' ...
         };
scale04={'off_nad_ang_platform'};
scale05={};
scale06={'swh_sq','hi_swh_sq'};
scale07={'lat','hi_lat','lon','hi_lon'};
scale09={'swh_sq_sd'};
map20={'flag_range_av','flag_range_ocog_av','flag_swh_sq_av',...
       'flag_sigma0_av','flag_sigma0_ocog_av','retrk_flag' ...
       };
map40={};
map80={};

%% translate all variables to scaled versions, checking for default values
fldnms=fieldnames(v);

for i=1:length(fldnms)
    fldnm=char(fldnms(i));
    switch fldnm
        case noscale
            v.(fldnm) = v_in.(fldnm);
        case scale01
            v.(fldnm) = double(v_in.(fldnm)).*1e-1;
        case scale02
            v.(fldnm) = double(v_in.(fldnm)).*1e-2;
        case scale03
            v.(fldnm) = double(v_in.(fldnm)).*1e-3;
        case scale04
            v.(fldnm) = double(v_in.(fldnm)).*1e-4;
        case scale05
            v.(fldnm) = double(v_in.(fldnm)).*1e-5;
        case scale06
            v.(fldnm) = double(v_in.(fldnm)).*1e-6;
        case scale07
            v.(fldnm) = double(v_in.(fldnm)).*1e-7;
        case scale09
            v.(fldnm) = double(v_in.(fldnm)).*1e-9;
        case map20
            if strcmp(fldnm,'retrk_flag')
                v.(fldnm) = map20_flag_unpack(v_in.(fldnm));
            else
                v.(fldnm) = map20_flag_fdm_unpack(v_in.(fldnm));
            end
        case map40
            v.(fldnm) = map40_flag_unpack(v_in.(fldnm));
        case map80
            v.(fldnm) = map80_flag_unpack(v_in.(fldnm));
        case 'flag_mcd'
            v.(fldnm).block_degraded = bitget(v_in.meas_conf_wd,32-0);
            v.(fldnm).blank_block = bitget(v_in.meas_conf_wd,32-1);
            v.(fldnm).datation_degrad = bitget(v_in.meas_conf_wd,32-2);
            v.(fldnm).orbit_propag_err = bitget(v_in.meas_conf_wd,32-3);
            v.(fldnm).orbit_file_change = bitget(v_in.meas_conf_wd,32-4);
            v.(fldnm).orbit_discontinuity = bitget(v_in.meas_conf_wd,32-5);
            v.(fldnm).echo_saturation = bitget(v_in.meas_conf_wd,32-6);
            v.(fldnm).other_echo_err = bitget(v_in.meas_conf_wd,32-7);
            v.(fldnm).rx_ch1_err = bitget(v_in.meas_conf_wd,32-8);
            v.(fldnm).rx_ch2_err = bitget(v_in.meas_conf_wd,32-9);
            v.(fldnm).window_delay_incon = bitget(v_in.meas_conf_wd,32-10);
            v.(fldnm).agc_incon = bitget(v_in.meas_conf_wd,32-11);
            v.(fldnm).cal1_corr_miss = bitget(v_in.meas_conf_wd,32-12);
            v.(fldnm).cal1_corr_ipf = bitget(v_in.meas_conf_wd,32-13);
            v.(fldnm).doris_uso_corr = bitget(v_in.meas_conf_wd,32-14);
            v.(fldnm).complex_cal1_ipf = bitget(v_in.meas_conf_wd,32-15);
            v.(fldnm).trk_echo_err = bitget(v_in.meas_conf_wd,32-16);
            v.(fldnm).rx1_echo_err = bitget(v_in.meas_conf_wd,32-17);
            v.(fldnm).rx2_echo_err = bitget(v_in.meas_conf_wd,32-18);
            v.(fldnm).noise_pow_incon = bitget(v_in.meas_conf_wd,32-19);
            v.(fldnm).azim_cal_miss = bitget(v_in.meas_conf_wd,32-20);
            v.(fldnm).default_azim_cal_ipf = bitget(v_in.meas_conf_wd,32-21);
            v.(fldnm).range_cal_miss = bitget(v_in.meas_conf_wd,32-22);
            v.(fldnm).range_cal_ipf = bitget(v_in.meas_conf_wd,32-23);
            v.(fldnm).cal2_corr_miss = bitget(v_in.meas_conf_wd,32-25);
            v.(fldnm).cal2_corr_ipf = bitget(v_in.meas_conf_wd,32-26);
            v.(fldnm).power_scaling_err = bitget(v_in.meas_conf_wd,32-27);
            v.(fldnm).attitude_corr_miss = bitget(v_in.meas_conf_wd,32-28);
            v.(fldnm).attitude_corr_err = bitget(v_in.meas_conf_wd,32-29);
    end
    switch fldnm
        case [notflag, scale01, scale02, scale03, scale04, scale05, scale06, scale07, scale09]
            %fprintf(['Field ' fldnm ' class ' class(v_in.(fldnm)) ' intmax = %d\n'],intmax(class(v_in.(fldnm))));
            v.(fldnm)(v_in.(fldnm)>=intmax(class(v_in.(fldnm)))) = NaN;
    end
end