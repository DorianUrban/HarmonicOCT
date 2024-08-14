function scanDetails = readLogFile(filename, folder)
%READLOGFILE Read the log file from a recording created by the Research
%Capture Software
scanDetails = struct;
scanDetails.is_c_scan = false;
scanDetails.ascans_per_cscan = 0;

opts = delimitedTextImportOptions;
opts.Delimiter = ': ';

filename = strcat(filename, folder, '\', folder, ' scan details.log');

try
    t = readtable(filename, opts);
catch
    error(strcat('Log file not found.', '\n', filename))
end

scanDetails.has_back_scan = false;

if sum(t.Var1 == "Scan type") == 1
    scanDetails.scan_shape = t.ExtraVar1{t.Var1 == "Scan type"};
    ymd = filename{1}(end-34:end-27);
    if strcmpi(scanDetails.scan_shape, "Line")
        if ymd>=20200709
            scanDetails.has_back_scan = false;
        else
            scanDetails.has_back_scan = true;
        end
    end
    if strcmpi(scanDetails.scan_shape, "Circular")
        scanDetails.has_back_scan = false;
    end
end

if sum(t.Var1 == "A-scan length") == 1
    scanDetails.ascan_length = str2num(t.ExtraVar1{t.Var1 == "A-scan length"}); %#ok<ST2NM>
else
    error('A-scan length not found')
end
if sum(t.Var1 == "Repeats per elevation") == 1
    scanDetails.raster_repeats = str2num(t.ExtraVar1{t.Var1 == "Repeats per elevation"}); %#ok<ST2NM>
else
    scanDetails.raster_repeats = 1;
end

ascan_per_scan_found = false;
if sum(t.Var1 == "C-scan dimensions (in A-scans)") == 1
    temp = t.ExtraVar1{t.Var1 == "C-scan dimensions (in A-scans)"};
    cscan_dimensions = split(temp, " high by ");
    scanDetails.cscan_h = str2num(cscan_dimensions{2}(1:end - 5));
    scanDetails.cscan_v = str2num(cscan_dimensions{1});
    scanDetails.is_c_scan = true;
elseif sum(t.Var1 == "Scan dimensions (in A-scans)") == 1
    temp = t.ExtraVar1{t.Var1 == "Scan dimensions (in A-scans)"};
    OCTA_dimensions = split(temp, " high by ");
    scanDetails.cscan_h = str2num(OCTA_dimensions{2}(1:end - 5));
    scanDetails.cscan_v = str2num(OCTA_dimensions{1});
    scanDetails.ascans_per_scan = (scanDetails.cscan_h * scanDetails.cscan_v * scanDetails.raster_repeats) + 200;
    ascan_per_scan_found = true;
    scanDetails.has_back_scan = true;
    scanDetails.return_path_ascans = 200;
elseif sum(t.Var1 == "OCTA raster dimensions (in A-scans)") == 1
    temp = t.ExtraVar1{t.Var1 == "OCTA raster dimensions (in A-scans)"};
    OCTA_dimensions = split(temp, " high by ");
    scanDetails.cscan_h = str2num(OCTA_dimensions{2}(1:end - 5));
    scanDetails.cscan_v = str2num(OCTA_dimensions{1});
    scanDetails.ascans_per_scan = (scanDetails.cscan_v * scanDetails.cscan_h * scanDetails.raster_repeats) + 200;
    ascan_per_scan_found = true;
    scanDetails.has_back_scan = true;
    scanDetails.return_path_ascans = 200;
end

% if sum(t.Var1 == "B-scan length (in A-scans)") == 1
%     scanDetails.ascans_per_scan = round(str2num(t.ExtraVar1{t.Var1 == "B-scan length (in A-scans)"})); %#ok<ST2NM>
%     ascan_per_scan_found = true;
% end

if sum(t.Var1 == "A-scans per B-scan") == 1
    scanDetails.ascans_per_scan = round(str2num(t.ExtraVar1{t.Var1 == "A-scans per B-scan"})); %#ok<ST2NM>
    ascan_per_scan_found = true;
end

if sum(t.Var1 == "A-scans per scan") == 1
    scanDetails.ascans_per_scan = round(str2num(t.ExtraVar1{t.Var1 == "A-scans per scan"})); %#ok<ST2NM>
    ascan_per_scan_found = true;
end

if sum(t.Var1 == "A-scans per C-scan") == 1
    scanDetails.ascans_per_scan = round(str2num(t.ExtraVar1{t.Var1 == "A-scans per C-scan"})); %#ok<ST2NM>
    ascan_per_scan_found = true;
end

if sum(t.Var1 == "Number of a-scans in back scan") == 1
    scanDetails.return_path_ascans = str2num(t.ExtraVar1{t.Var1 == "Number of a-scans in back scan"});
else
    if sum(t.Var1 == "B-scan length (in A-scans)")
        scanDetails.return_path_ascans = scanDetails.ascans_per_scan - str2num(t.ExtraVar1{t.Var1 == "B-scan length (in A-scans)"});
    else
        scanDetails.return_path_ascans = 0;
    end
end

%     error('Could not find either a-scans per b-scan or a-scans per c-scan.')
% end

if sum(t.Var1 == "Scan duration") == 1
    scan_duration = t.ExtraVar1{t.Var1 == "Scan duration"};
    scanDetails.scan_duration = str2double(scan_duration);
else
    error('Scan duration not found.')
end


if sum(t.Var1 == "Analogue output (galvo) sample rate") == 1
    AO_sample_rate  = t.ExtraVar1{t.Var1 == "Analogue output (galvo) sample rate"};
    AO_sample_rate = AO_sample_rate(1:end-2);
    scanDetails.AO_sample_rate = str2double(AO_sample_rate);
else
    error('Analogue input sample rate not found.')
end

if sum(t.Var1 == "Analogue input (galvo and APD) sample rate") == 1
    AI_sample_rate  = t.ExtraVar1{t.Var1 == "Analogue input (galvo and APD) sample rate"};
    AI_sample_rate = AI_sample_rate(1:end-2);
    scanDetails.AI_sample_rate = str2double(AI_sample_rate);
else
    if sum(t.Var1 == "Analogue input (galvo) sample rate") == 1
        AI_sample_rate  = t.ExtraVar1{t.Var1 == "Analogue input (galvo) sample rate"};
        AI_sample_rate = AI_sample_rate(1:end-2);
        scanDetails.AI_sample_rate = str2double(AI_sample_rate);
    else
        error("Analogue input galvo sample rate not found.")
    end        
end

if sum(t.Var1 == "Number of start park scans") == 1
    scanDetails.pre_scans = str2num(t.ExtraVar1{t.Var1 == "Number of start park scans"});
end

if sum(t.Var1 == "Number of end park scans") == 1
    scanDetails.post_scans = str2num(t.ExtraVar1{t.Var1 == "Number of end park scans"});
else
    scanDetails.post_scans = 0;
end

if sum(t.Var1 == "Number of start park scans") == 1
    scanDetails.pre_scans = str2num(t.ExtraVar1{t.Var1 == "Number of start park scans"}); %#ok<ST2NM>
end

if sum(t.Var1 == "Number of end park scans") == 1
    scanDetails.post_scans = str2num(t.ExtraVar1{t.Var1 == "Number of end park scans"}); %#ok<ST2NM>
end

scanDetails.FFT_output_length = scanDetails.ascan_length/2;
if sum(t.Var1 == "Number of scans captured") == 1
    scanDetails.number_of_scans = str2num(t.ExtraVar1{t.Var1 == "Number of scans captured"}); %#ok<ST2NM>
    
else
%     warning('Number of frames not found.')
%     error('Number of frames not found.')
    scanDetails.number_of_scans = 1;
end


if sum(t.Var1 == "Number of return path a-scans")
    scanDetails.flyback_ascans = str2num(t.ExtraVar1{t.Var1 == "Number of return path a-scans"}); %#ok<ST2NM>
end

if sum(t.Var1 == "APD data") == 1
    test = t.ExtraVar1(t.Var1 == "APD data");
    if strcmp(test, "Yes") == true
        scanDetails.APD_data = true;
    else
        scanDetails.APD_data = false;
    end
else
    scanDetails.APD_data = true;
end

if sum(t.Var1 == "OCT source type") == 1
    scanDetails.source_type = t.ExtraVar1(t.Var1 == "OCT source type");
    scanDetails.source_type = scanDetails.source_type{1};
	else
	scanDetails.source_type = "unknown";
end

if sum(t.Var1 == "Lateral resolution (micrometres)") == 1
    scanDetails.lateralResolution = t.ExtraVar1(t.Var1 == "Lateral resolution (micrometres)");
    scanDetails.lateralResolution = str2num(scanDetails.lateralResolution{1}) * 1e-6;
end

if sum(t.Var1 == "Dispersion compensation setting") == 1
    scanDetails.dispersionCompensation = t.ExtraVar1(t.Var1 == "Dispersion compensation setting");
    scanDetails.dispersionCompensation = str2num(scanDetails.dispersionCompensation{1});
end


if sum(t.Var1 == "A-scan sample rate") == 1
    temp = t.ExtraVar1(t.Var1 == "A-scan sample rate");
    temp = temp{1}(1:end-2);
    temp = str2num(temp);
    scanDetails.ascanRate = temp;
end

scanDetails.comment = t.Var1{5};

if sum(t.Var1 == "Spectrometer coefficients") == 1
    temp = t.ExtraVar1(t.Var1 == "Spectrometer coefficients");
    temp = split(temp, ', ');
    temp = temp(1:end-1);
    for count=1:numel(temp)
        temp{count} = str2num(temp{count});
    end
    temp = cell2mat(temp);
    scanDetails.spectrometer_coefficients = temp;
end


% We need flyback_ascans for a-scans

if ~isfield(scanDetails,'flyback_ascans')
    if strcmpi(scanDetails.scan_shape, "H-scan with scan pattern phase shift")
        scanDetails.flyback_ascans = scanDetails.ascans_per_scan/2;
    elseif strcmpi(scanDetails.scan_shape, "V-scan with scan pattern phase shift")
        scanDetails.flyback_ascans = scanDetails.ascans_per_scan/2;
    else
        scanDetails.flyback_ascans = 0;
    end
end


if sum(t.Var1 == "Stimulus mode") == 1
    temp = t.ExtraVar1(t.Var1 == "Stimulus mode");   
    scanDetails.stimulus_mode = temp;
scanDetails.stimulus_mode = (t.ExtraVar1{t.Var1 == "Stimulus mode"});

end

if sum(t.Var1 == "Stimulus delay (scans)") == 1
    scanDetails.stimulus_delay = str2num(t.ExtraVar1{t.Var1 == "Stimulus delay (scans)"});
end
 
if sum(t.Var1 == "Stimulus duration (ms)") == 1
    scanDetails.stimulus_duration = str2num(t.ExtraVar1{t.Var1 == "Stimulus duration (ms)"});
end

%%% Dorian Added

if any(t.Var1 == "Focus motor position")
    scanDetails.focus_position = str2num(t.ExtraVar1{t.Var1 == "Focus motor position"});
end

if any(t.Var1 == "Reference motor position")
    scanDetails.reference_position = str2num(t.ExtraVar1{t.Var1 == "Reference motor position"});
end

end

