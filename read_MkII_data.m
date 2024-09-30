function scan_data = read_MkII_data(path)
    persistent scan_data_path;
    fid = fopen(path, 'r');
    
    if isequal(fid, -1)
        [path,name,ext] = fileparts(path);
        if isempty(path) && (isempty(ext) || strcmp(ext,'.txt'))
            if isempty(scan_data_path)
                scan_data_path = fullfile(fileparts(mfilename('fullpath')), ...
                    '../../../SHeM Data/Scan Experiments/Scan Data/');
            end
            scan_data = read_MkII_data(fullfile(scan_data_path, [name '.txt']));
            return
        else
            error(['Failed to open: ' path])
        end
    end
    
    scan_data = struct();
    
    str = fgetl(fid);
    scan_data.start_time = datetime(str(28:end), 'InputFormat', 'd/MM/yyyy h:mm:ss a');
    
    str = fgetl(fid);
    scan_data.end_time = datetime(str(26:end), 'InputFormat', 'd/MM/yyyy h:mm:ss a');
    
    % Start position
    start = [0 0 0];
    start(3) = fscanf(fid, '%%Sample Z (um): %f'); fseek(fid,2,0);
    start(1) = fscanf(fid, '%%Initial X (um): %f'); fseek(fid,2,0);
    start(2) = fscanf(fid, '%%Initial Y (um): %f'); fseek(fid,2,0);
    scan_data.start = start;
    
    % Scan size
    sz = [0 0];
    sz(1) = fscanf(fid, '%%Scan size X (um): %f'); fseek(fid,2,0);
    sz(2) = fscanf(fid, '%%Scan size Y (um): %f'); fseek(fid,2,0);
    scan_data.image_size = sz;
    
    % Number of pixels
    pixels = [0 0];
    pixels(1) = fscanf(fid, '%%Number of X pixels: %d'); fseek(fid,2,0);
    pixels(2) = fscanf(fid, '%%Number of Y pixels: %d'); fseek(fid,2,0);
    scan_data.num_pixels = pixels;
    
    scan_data.step_size = sz ./ pixels;
    
    for i = 1:4
        fgetl(fid);
    end
    
    scan_data.read_time = fscanf(fid, '%%Detector read time (ms): %f'); fseek(fid,2,0);
    scan_data.wait_time = fscanf(fid, '%%Detector wait time (ms): %f'); fseek(fid,2,0);
    
    fgetl(fid);
    fseek(fid, 16, 0);
    
    scan_data.comments = fgetl(fid);
    
    scan_data.image = fscanf(fid, '%f', pixels);
    
%     try
%         scan_data.system = get_system_log_period(scan_data.start_time, scan_data.end_time);
%     catch
%         warning('Failed to read system logs, potentially using old file format')
%     end
    
    fclose(fid);
end