function [frame, scanDetails, ch2_data] = get_frame(folder, number, mean_subtract)
    
%     if nargin < 5, ch2 = false; end
%     if nargin < 4, complex = false; end
    if nargin < 3, mean_subtract = false; end
    if nargin < 2, number = 1; end

    path = split(folder, "\");
    file = string(path(end));
    directory = strjoin(path(1:end - 1), "\");
    directory = strcat(directory, "\");

    scanDetails = readLogFile(directory, file);

    all_files = dir(folder);
    all_files = {all_files(3:end).name};


    data = regexp(all_files, "\d+(?=.bin)", "match");


    % % start_num = min(str2num(cell2mat([data{:}].')));
    % start_num = min(cellfun(@str2num, [data{:}].'));
    % number = number + start_num;

    idx = ~cellfun("isempty", regexp(all_files, "Alazartech OCT", "match"));
    data_files = all_files(idx);
    % filepath = strcat(folder, "\", file, sprintf("- Alazartech OCT %d.bin", number));
    filepath = strcat(folder, "\", data_files{number});

    file = fopen(filepath, 'r');
    data = fread(file, [2 * scanDetails.ascan_length Inf], 'uint16');
    fclose(file);
    
    frame = data(1:2:end, :);

    if mean_subtract
        frame = frame - mean(frame, 2);
        frame = frame - movmean(frame, 201);
    end
    
    if nargout > 2
        ch2_data = data(2:2:end, :);
    end
end