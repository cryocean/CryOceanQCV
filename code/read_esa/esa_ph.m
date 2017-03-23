function [HDR] = esa_ph(PH)

% Function to read an Product Header (PH) from an
% ESA data File into a structure only retaining useful records
% within the MPH, ignoring blanks, quotes and newlines
% MPH is a character string with entries such as:
%       STATE_VECTOR_TIME="01-JUL-2013 03:07:44.940476"
%       DELTA_UT1=+.000000<s>
% With newline delimeters
% Structure item names taken from entries
% Where units are provided, an additional field 'title'_units is generated
% DataSet records will be placed in a DS structure within HDR
%   HDR.DS(1:n) where n is the number of dataset records included
% [HDR] = esa_ph(PH)
% Valid for ENVISAT and Cryosat-2 Products

%% Set dataet flag off and number of dataset read = 0
is_ds = 0;
k=0;

%% Split product header into lines
PHL = regexp(PH,'\n','split');

%% Save data to Structure

for i=1:length(PHL)
    
    tline = PHL{i};
    
    % equal sign delimits end of header field name
    I=strfind(tline,'=')+1;
    fld=lower(strtrim(tline(1:I-2)));
    
    if ~isempty(fld)
        % For DataSet entries - should put in HDR.DS(k).
        if strcmp(fld,'ds_name')
            is_ds = 1;
            k=k+1;
        end
        
        % For pure text entry - header entry delimited by ""
        if strcmp(tline(I),'"')
            value = strtrim(tline(I+1:end-1));
            units = '';
            scale = 1.;
            
        else
            % Units will be delimited by <>
            J=strfind(tline,'<')-1;
            
            % If no units - use entire string
            if isempty(J)
                J=length(tline);
                units = '';
                scale = 1.;
            else
                % Otherwise, save units
                J2 = strfind(tline,'>')-1;
                units = tline(J+2:J2);
            end
            
            if not(isempty(tline(I:J)))
                if isnan(str2double(tline(I:J)))
                    value=strtrim(tline(I:J));
                else
                    value=str2double(tline(I:J));
                end
            end
        end
        
        if ~isempty(value)
            if is_ds
                HDR.DS(k).(fld) = value;
                if ~isempty(units)
                    HDR.DS(k).([fld '_units']) = units;
                end
            else
                HDR.(fld) = value;
                if ~isempty(units)
                    HDR.([fld '_units']) = units;
                end
            end
        end
    end
end
return;

