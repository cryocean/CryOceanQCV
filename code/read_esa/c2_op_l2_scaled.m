function [v] = c2_op_l2_scaled(v_in,n_recs)

% Scale the SIR_L2_IOP / SIR_L2_GOP records from a
%  CRYOSAT L2 (I/G)OP File, stored by variable,
%  to SI units, setting default values to NaNs
% 'map' flags converted to arrays
% [v] = c2_op_l2_scaled(v_in,n_recs)
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
    ...%'tai_days', cast(zeros(n_recs,1), 'int32'), ...      % TAI Time stamp of MDSR - days
    ...%'tai_secs', cast(zeros(n_recs,1), 'uint32'), ...     % TAI Time stamp of MDSR - seconds
    ...%'tai_usecs', cast(zeros(n_recs,1), 'uint32'), ...    % TAI Time stamp of MDSR - 10-6 secs
    'time_offset', cast(zeros(n_recs,1), 'int16'), ...       % TAI - UTC offset
    'utc_time', zeros(n_recs,1), ...                         % full UTC time
    'tdiff', zeros(n_recs,n_hi), ...          % hi rate time difference of hi rate from 1Hz [n_hi]
    'hi_time', zeros(n_recs,n_hi), ...                       % hi rate full TAI [n_hi]
    'tdiff_offset', cast(zeros(n_recs,n_hi), 'int16'), ...   % hi rate TAI - UTC offset [n_hi]
    'hi_utctime', zeros(n_recs,n_hi), ...                    % hi rate full UTC time [n_hi]
    'rec_cnt', cast(zeros(n_recs,1), 'uint32'), ...          % Record Counter
    'lat', zeros(n_recs,1), ...               % Geodetic Latitude (positive N)
    'hi_lat', zeros(n_recs,n_hi), ...         % hi rate Geodetic Latitude (positive N) [n_hi]
    'lon', zeros(n_recs,1), ...               % Longitude (positive E, 0 at Greenwich)
    'hi_lon', zeros(n_recs,n_hi), ...         % hi rate Longitude (positive E, 0 at Greenwich) [n_hi]
    'h', zeros(n_recs,1), ...                  % Altitude of CoG above reference ellipsoid
    'hi_h', zeros(n_recs,n_hi), ...            % hi rate altitude differences from 1 Hz altitude [n_hi].
    'h_rate', zeros(n_recs,1), ...             % Instantaneous altitude rate
    'meas_conf_wd', cast(zeros(n_recs,n_hi),'uint32'), ...      % Measurement Confidence Word
    'flag_mcd', struct( ...
          'block_degraded', cast(zeros(n_recs,n_hi),'int8'), ...
          'blank_block', cast(zeros(n_recs,n_hi),'int8'), ...
          'orbit_propag_err', cast(zeros(n_recs,n_hi),'int8'), ...
          'orbit_file_change', cast(zeros(n_recs,n_hi),'int8'), ...
          'orbit_discontinuity', cast(zeros(n_recs,n_hi),'int8'), ...
          'echo_saturation', cast(zeros(n_recs,n_hi),'int8'), ...
          'other_echo_err', cast(zeros(n_recs,n_hi),'int8'), ...
          'cal1_corr_miss', cast(zeros(n_recs,n_hi),'int8'), ...
          'cal1_corr_ipf', cast(zeros(n_recs,n_hi),'int8'), ...
          'doris_uso_corr', cast(zeros(n_recs,n_hi),'int8'), ...
          'trk_echo_err', cast(zeros(n_recs,n_hi),'int8'), ...
          'rx1_echo_err', cast(zeros(n_recs,n_hi),'int8'), ...
          'rx2_echo_err', cast(zeros(n_recs,n_hi),'int8'), ...
          'cal2_corr_miss', cast(zeros(n_recs,n_hi),'int8'), ...
          'cal2_corr_ipf', cast(zeros(n_recs,n_hi),'int8'), ...
          'power_scaling_err', cast(zeros(n_recs,n_hi),'int8'), ...
          'processing_type', cast(zeros(n_recs,n_hi),'int8') ...
         ), ...
    'pp', zeros(n_recs,1), ...                 % Pulse Peakiness
    'hi_pp', zeros(n_recs,n_hi), ...           % hi rate Pulse Peakiness [n_hi]
    'trk_error', zeros(n_recs,n_hi), ...       % hi rate tracker error [n_hi]
    'flag_trk', cast(zeros(n_recs,n_hi),'uint8'), ...       % hi rate tracker flag [n_hi]
    ...% Range Measurements Group
    'range', zeros(n_recs,1), ...             % ocean range
    'hi_range', zeros(n_recs,n_hi), ...       % hi rate ocean range [n_hi]
    'range_sd', zeros(n_recs,1), ...          % Standard deviation of hi rate ocean range
    'range_num', cast(zeros(n_recs,1),'uint16'), ...         % Number of valid hi rate points for ocean range
    'flag_range_av', cast(zeros(n_recs,n_hi),'uint8'), ...     % ocean range averaging flag
    'range_ice', zeros(n_recs,1), ...         % ice ranges
    'hi_range_ice', zeros(n_recs,n_hi), ...   % hi rate ice2 ranges [n_hi]
    'range_ice_sd', zeros(n_recs,1), ...      % Standard deviation of hi rate ice ranges
    'range_ice_num', cast(zeros(n_recs,1),'uint16'), ...     % Number of valid hi rate points for ice ranges
    'flag_range_ice_av', cast(zeros(n_recs,n_hi),'uint8'), ... % ice range averaging flag
    ...% Range Corrections Group
    'corr_dopp', cast(zeros(n_recs,1),'int16'), ...          % Doppler correction
    'corr_uso', cast(zeros(n_recs,1),'int16'), ...           % range correction for USO drift
    'corr_cog', cast(zeros(n_recs,1),'int16'), ...        % range correction for Centre of Gravity
    'corr_int_cal1', cast(zeros(n_recs,1),'int16'), ...      % range correction for internal CAL1
    'corr_instr', cast(zeros(n_recs,1),'int16'), ...	     % range instrumental correction
    'corr_dry_trop', cast(zeros(n_recs,1),'int16'), ...	     % Model dry tropospheric correction
    'corr_wet_trop_mod', cast(zeros(n_recs,1),'int16'), ...  % Model wet tropospheric correction
    'corr_inv_baro', cast(zeros(n_recs,1),'int16'), ...      % Inverted barometer correction
    'corr_dac', cast(zeros(n_recs,1),'int16'), ...           % Dynamic Atmospheric Correction
    'corr_iono_gim', cast(zeros(n_recs,1),'int16'), ...      % Ionospheric correction from GIM model
    'corr_ssb', cast(zeros(n_recs,1),'int16'), ...           % Sea state bias
    ...% SWH & Backscatter Measurements Group
    'swh_sq', zeros(n_recs,1), ...             % square of swh
    'swh', zeros(n_recs,1), ...                % Significant wave height
    'hi_swh', zeros(n_recs,n_hi), ...          % hi rate swh [n_hi]
    'swh_sd', zeros(n_recs,1), ...            % Standard deviation of hi rate swh
    'swh_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for swh
    'flag_swh_av', cast(zeros(n_recs,n_hi),'uint8'), ...       % swh averaging flag
    'sigma0', zeros(n_recs,1), ...                % corrected Ocean backscatter coefficient
    'hi_sigma0', zeros(n_recs,n_hi), ...          % hi rate sigma0 [n_hi]
    'sigma0_sd', zeros(n_recs,1), ...            % Standard deviation of hi rate sigma0
    'sigma0_num', cast(zeros(n_recs,1),'uint16'), ...           % Number of valid hi rate points for sigma0
    'flag_sigma0_av', cast(zeros(n_recs,n_hi),'uint8'), ...       % sigma0 averaging flag
    'sigma0_ice', zeros(n_recs,1), ...            % corrected Ice backscatter coefficient
    'hi_sigma0_ice', zeros(n_recs,n_hi), ...      % hi rate ice sigma0 [n_hi]
    'sigma0_ice_sd', zeros(n_recs,1), ...        % Standard deviation of hi rate ice sigma0
    'sigma0_ice_num', cast(zeros(n_recs,1),'uint16'), ...       % Number of valid hi rate points for ice sigma0
    'flag_sigma0_ice_av', cast(zeros(n_recs,n_hi),'uint8'), ...   % ice sigma0 averaging flag
    'off_nad_ang_sq_wvform', zeros(n_recs,1), ... % Off nadir angle squared of the satellite from waveform data
    'agc', zeros(n_recs,1), ...                   % AGC
    'hi_sigma0_scale',zeros(n_recs,n_hi), ...    % hi rate scale factor for sigma0 [n_hi]
    'corr_swh_instr', zeros(n_recs,1), ...        % SWH instrument correction
    'corr_agc', zeros(n_recs,1), ...              % AGC correction
    'corr_sigma0_int_cal1', zeros(n_recs,1), ...  % internal cal1 correction on sigma0
    'corr_sigma0_instr', zeros(n_recs,1), ...     % instrument correction on sigma0
    'corr_atm_atten', zeros(n_recs,1), ...        % atmospheric attenuation correction on sigma0
    ...% Geophysical Corrections Group
    'h_mss1', zeros(n_recs,1), ...                % Mean sea-surface height 1
    'h_mss2', zeros(n_recs,1), ...                % Mean sea-surface height 2
    'h_geoid', zeros(n_recs,1), ...               % Geoid height
    'h_dem', zeros(n_recs,1), ...                 % Ocean depth/land elevation (DEM)
    'h_mdt', zeros(n_recs,1), ...                 % Mean dynamic topography
    'h_tot_geocen_ocn_tide_sol1', zeros(n_recs,1), ...	% Total geocentric ocean tide height (solution 1)
    'h_tot_geocen_ocn_tide_sol2', zeros(n_recs,1), ...	% Total geocentric ocean tide height (solution 2)
    'h_long_period_ocn_tide', zeros(n_recs,1), ...	    % Long period Tide height
    'h_long_period_ocn_tide_no_eq', zeros(n_recs,1), ...	% Long period Tide height - none equlibrium tides
    'h_tide_load_sol1', zeros(n_recs,1), ...         % Tidal loading height (solution 1)
    'h_tide_load_sol2', zeros(n_recs,1), ...         % Tidal loading height (solution 2)
    'h_tide_solid', zeros(n_recs,1), ...             % Solid earth tide height
    'h_tide_geocen_pole', zeros(n_recs,1), ...       % Geocentric pole tide height
    'wsp', zeros(n_recs,1), ...               % RA2 wind speed 
    'wsp_u_mod', zeros(n_recs,1), ...         % u component of the model wind vector
    'wsp_v_mod', zeros(n_recs,1), ...         % v component of the model wind vector 
    'flag_surface_type', cast(zeros(n_recs,1),'uint16') ... % Surface type flag
);

%% Define which variables are left as they are, scaled or mapped
noscale={'time','time_offset','utc_time',...
          'hi_time','tdiff_offset','hi_utctime', ...
          'rec_cnt','meas_conf_wd','range_num','range_ice_num', ...
          'swh_num','sigma0_num','sigma0_ice_num', 'flag_surface_type'...
         };
notflag={'time_offset','tdiff_offset','rec_cnt','range_num','range_ice_num' ...
         'swh_num','sigma0_num','sigma0_ice_num' ...
         };
scale01={};
scale02={'pp','hi_pp',...
         'sigma0','hi_sigma0','sigma0_sd', ...
         'sigma0_ice','hi_sigma0_ice','sigma0_ice_sd',...
         'agc','hi_sigma0_scale','corr_agc',...
         'corr_sigma0_int_cal1','corr_sigma0_instr','corr_atm_atten' ...
         };
scale03={'h','hi_h','h_rate',...
         'range','hi_range','range_sd', ...
         'range_ice','hi_range_ice','range_ice_sd', ...
         'corr_dopp','corr_uso','corr_cog','corr_int_cal1','corr_instr', ...
         'corr_dry_trop','corr_wet_trop_mod','corr_inv_baro','corr_dac', ...
         'corr_iono_gim','corr_ssb', ...
         'swh','hi_swh','swh_sd','corr_swh_instr', ...
         'h_mss1','h_mss2','h_geoid','h_dem','h_mdt', ...
         'h_tot_geocen_ocn_tide_sol1','h_tot_geocen_ocn_tide_sol2', ...
         'h_long_period_ocn_tide','h_long_period_ocn_tide_no_eq', ...
         'h_tide_load_sol1','h_tide_load_sol2',...
         'h_tide_solid','h_tide_geocen_pole',...
         'wsp','wsp_u_mod','wsp_v_mod' ...
         };
scale04={'trk_error','off_nad_ang_sq_wvform'};
scale05={};
scale06={'tdiff','swh_sq'};
scale07={'lat','hi_lat','lon','hi_lon'};
map20={'flag_trk','flag_range_av','flag_range_ice_av','flag_swh_av',...
       'flag_sigma0_av','flag_sigma0_ice_av'...
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
        case map20
            v.(fldnm) = map20_flag_unpack(v_in.(fldnm));
        case map40
            v.(fldnm) = map40_flag_unpack(v_in.(fldnm));
        case map80
            v.(fldnm) = map80_flag_unpack(v_in.(fldnm));
        case 'flag_mcd'
            v.(fldnm).block_degraded = bitget(v_in.meas_conf_wd,32-0);
            v.(fldnm).blank_block = bitget(v_in.meas_conf_wd,32-1);
            v.(fldnm).orbit_propag_err = bitget(v_in.meas_conf_wd,32-3);
            v.(fldnm).orbit_file_change = bitget(v_in.meas_conf_wd,32-4);
            v.(fldnm).orbit_discontinuity = bitget(v_in.meas_conf_wd,32-5);
            v.(fldnm).echo_saturation = bitget(v_in.meas_conf_wd,32-6);
            v.(fldnm).other_echo_err = bitget(v_in.meas_conf_wd,32-7);
            v.(fldnm).cal1_corr_miss = bitget(v_in.meas_conf_wd,32-12);
            v.(fldnm).cal1_corr_ipf = bitget(v_in.meas_conf_wd,32-13);
            v.(fldnm).doris_uso_corr = bitget(v_in.meas_conf_wd,32-14);
            v.(fldnm).trk_echo_err = bitget(v_in.meas_conf_wd,32-16);
            v.(fldnm).rx1_echo_err = bitget(v_in.meas_conf_wd,32-17);
            v.(fldnm).rx2_echo_err = bitget(v_in.meas_conf_wd,32-18);
            v.(fldnm).cal2_corr_miss = bitget(v_in.meas_conf_wd,32-25);
            v.(fldnm).cal2_corr_ipf = bitget(v_in.meas_conf_wd,32-26);
            v.(fldnm).power_scaling_err = bitget(v_in.meas_conf_wd,32-27);
            v.(fldnm).processing_type = ...
              bitget(v_in.meas_conf_wd,32-28) + 2.*bitget(v_in.meas_conf_wd,32-29);
    end
    switch fldnm
        case [notflag, scale01, scale02, scale03, scale04, scale05, scale06, scale07]
            %fprintf(['Field ' fldnm ' class ' class(v_in.(fldnm)) ' intmax = %d\n'],intmax(class(v_in.(fldnm))));
            v.(fldnm)(v_in.(fldnm)>=intmax(class(v_in.(fldnm)))) = NaN;
    end
end