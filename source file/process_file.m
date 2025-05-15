
function process_file(currentFile, file_folder, parentFolder, params, file_idx, total_files)
% PROCESS_FILE - Processes a single spike data file
%
% INPUTS:
%   currentFile - File information from dir()
%   file_folder - Path to the folder containing the file
%   parentFolder - Path to the output folder
%   params - Parameter structure
%   file_idx - Current file index
%   total_files - Total number of files

    % Extract file information
    fileName = currentFile.name;
    baseFileName = fileName(1:end-4); % Remove .spk extension
    
    % Create output folder for this file
    outputFolder = fullfile(parentFolder, baseFileName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Log processing
    fprintf('Processing file %d/%d: %s\n', file_idx, total_files, fileName);
    
    % Load spike data
    filePath = fullfile(file_folder, fileName);
    try
        allData = AxisFile(filePath).SpikeData.LoadData;
        SPK_dataset = AxisFile(filePath).SpikeData;
        fr = SPK_dataset.SamplingFrequency; % Hz
        fprintf('Sampling frequency: %d Hz\n', fr);
        
        % Get data dimensions
        [nwr, nwc, nec, ner] = size(allData);
        fprintf('Data dimensions: [%d, %d, %d, %d]\n', nwr, nwc, nec, ner);
    catch ME
        fprintf('Error loading file: %s\nMessage: %s\n', filePath, ME.message);
        return;
    end
    
    % Initialize measurement tables
    % Initialize analysis tables
    [T, T_electrode, T_parameters] = initialize_tables(params);
    % Measurements = {'Number of Spikes', 'Mean Firing Rate (Hz)', 'ISI in general',...
    %     'Number of Bursts',...
    %     'Burst Frequency (Hz)',...
    %     'Burst Duration -Avg(sec)','Burst Duration -Std(sec)', 'Burst Duration -Median(sec)', ...
    %     'Burst Percentage (%)',...
    %     'Normalized Burst Duration IQR',...
    %     'Number of Spikes per Burst -Avg', 'Number of Spikes per Burst -Std', 'Number of Spikes per Burst -Median',...
    %     'Mean ISI within Burst -Avg', 'Mean ISI within Burst -Std', 'Mean ISI within Burst -Median',...
    %     'ISI Coefficient of Variation',...
    %     'Inter-Burst Interval - Avg', 'Inter-Burst Interval - Std', 'Inter-Burst Interval - Median',...
    %     'Inter-Burst Interval Coefficient of Variation',...
    %     'Amplitude -Avg(mV)','Amplitude -Std(mV)',...
    %     'Amplitude - Avg Peak to Trough(mV)', 'Amplitude - Std Peak to Trough(mV)',...
    %     'max ISI (ms)','number of removed short ISIs'};
    % T = table(Measurements');
    % 
    % Measurements_electrode = {'Total number of spikes', 'number of excluded spikes',...
    %     'unsorted firing rates (Hz) ', 'Number of Units Detected'};
    % T_electrode = table(Measurements_electrode');
    % 
    % % Save parameter list
    % parameter_list = {'threshold_to_merge', 'refractoryT', ...
    %     'min_spikes_E', 'max_ISI_E', 'min_spikes_N', 'max_ISI_N',...
    %     'network_participation_threshold',...
    %     'cutoff_frequency', 'flag_threshold', 'std_cutoff'};
    % T_parameters = table(parameter_list');
    % 
    % Tleft_parameters = table([params.threshold_to_merge; params.refractoryT;...
    %     params.min_spikes_E; params.max_ISI_E; params.min_spikes_N; params.max_ISI_N;...
    %     params.network_participation_threshold; ...
    %     params.cutoff_frequency; params.flag_threshold; params.std_cutoff]);
    % Tleft_parameters.Properties.VariableNames(1) = {'parameter Value'};
    % T_parameters = [T_parameters, Tleft_parameters];
    writetable(T_parameters, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'parameter list');
    
    % Initialize data structures
    raster_raw = {};
    sorting_results = {};
    num_well = 1;
    total_num_cell = 1;
    
    % Get recording properties
    [maxTime, maxSpikeAmp, minSpikeAmp] = get_recording_properties(allData);
    fprintf('Recording duration: %.2f seconds\n', maxTime);
    
    % Initialize statistics
    num_total_electrodes_with_spikes = 0;
    num_inactive_units = 0;
    num_total_detected_units = 0;
    fr_inactive_units = [];
    fr_active_units = [];
    
    % Initialize PowerPoint
    try
        pptx = exportToPPTX('SAMS1.pptx');
    catch ME
        warning('Error initializing PowerPoint: %s', ME.message);
        pptx = [];
    end
    
    % Initialize lists for potentially problematic electrodes
    possibleMUAList = {};
    flagTooMuchList = {};
    
    % Process each electrode
    for i = 1:nwr
        for j = 1:nwc
            num_electrode = 1;
            for m = 1:nec
                for n = 1:ner
                    % Update progress
                    frac4 = n/ner;
                    frac3 = ((m-1)+frac4)/nec;
                    frac2 = ((j-1)+frac3)/nwc;
                    frac1 = ((i-1)+frac2)/nwr;
                    frac0 = file_idx/total_files;
                    try
                        progressbar(frac0, frac1);
                    catch
                        % Continue without progress update
                    end
                    
                    % Get current electrode data
                    Spikes = allData{i, j, m, n}(:);
                    
                    % Skip if no spikes
                    if isempty(Spikes)
                        continue;
                    end
                    
                    % Get spike times and waveforms
                    [Times, spikes] = Spikes.GetTimeVoltageVector;
                    fprintf('Processing electrode %s: %d spikes\n', ...
                        [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)], size(spikes, 2));
                    
                    % Perform PCA
                    try
                        [~, score, ~] = pca(spikes');
                    catch ME
                        fprintf('Error in PCA: %s\n', ME.message);
                        continue;
                    end
                    
                    % Process electrode
                    try
                        [electrode_results, electrode_stats] = process_electrode(spikes, Times, i, j, m, n, ...
                            params, maxTime, maxSpikeAmp, minSpikeAmp, score);
                    catch ME
                        fprintf('Error processing electrode: %s\n', ME.message);
                        continue;
                    end
                    
                    % Skip if no units detected
                    if isempty(electrode_results.new_idx_list)
                        continue;
                    end
                    
                    % Update results
                    try
                        [T, T_electrode, raster_raw, sorting_results, num_electrode, total_num_cell, ...
                        num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                        fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList] = ...
                        update_results(electrode_results, electrode_stats, T, T_electrode, raster_raw, ...
                        sorting_results, num_electrode, total_num_cell, num_well, i, j, m, n, ...
                        num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                        fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList, pptx, ...
                        maxTime, maxSpikeAmp, minSpikeAmp, baseFileName, spikes, Times, score);
                    catch ME
                        fprintf('Error updating results: %s\n', ME.message);
                    end
                end
            end
            num_well = num_well + 1;
        end
    end
    % Update progress bar with completion message for this file
    try
        % Custom progress bar with completion message
        h = figure('Position', [400, 400, 400, 100], 'Name', 'Processing Status', ...
                  'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');
        uicontrol('Style', 'text', 'Position', [20, 50, 360, 30], ...
                 'String', ['File ' num2str(file_idx) '/' num2str(total_files) ' completed. Writing outputs...'], ...
                 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        uicontrol('Style', 'text', 'Position', [20, 20, 360, 20], ...
                 'String', 'Please wait. Do not close MATLAB.', ...
                 'FontSize', 10, 'HorizontalAlignment', 'center');
        drawnow;
    catch
        % If figure creation fails, just print a message
        fprintf('\nFile %d/%d completed. Writing outputs...\n', file_idx, total_files);
        fprintf('Please wait. Do not close MATLAB.\n');
    end
    % Save PowerPoint
    if ~isempty(pptx)
        try
            pptx.save([outputFolder, '\', baseFileName]);
            fprintf('PowerPoint presentation saved successfully.\n');
        catch ME
            warning('Error saving PowerPoint: %s', ME.message);
        end
    end
    
    % Save burst info
    save([outputFolder, '\burst_info_all.mat'], 'raster_raw', 'maxTime', 'sorting_results', '-v7.3');
    
    % Perform network burst analysis
    fprintf('Performing network burst analysis...\n');
    get_network_burst_info_08202024(raster_raw, maxTime, fr, params.network_participation_threshold, ...
        params.min_spikes_E, params.max_ISI_E, params.min_spikes_N, params.max_ISI_N, ...
        outputFolder, sorting_results);
    
    % Save results to Excel
    fprintf('Writing results to Excel...\n');
    writetable(T, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'individual unit');
    writetable(T_electrode, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'electrode statistics');

    % Save check lists
    if isempty(possibleMUAList)
        possibleMUAList{end+1} = 'no check list';
    end
    if isempty(flagTooMuchList)
        flagTooMuchList{end+1} = 'no check list';
    end

    % Make both tables have the same length
    maxLength = max(length(possibleMUAList), length(flagTooMuchList));
    while length(possibleMUAList) < maxLength
        possibleMUAList{end+1} = '';
    end
    while length(flagTooMuchList) < maxLength
        flagTooMuchList{end+1} = '';
    end

    % Create tables with equal numbers of rows
    T_pmua = cell2table(possibleMUAList', 'VariableNames', "Possible MultiUnit");
    T_flag = cell2table(flagTooMuchList', 'VariableNames', "over-exlcuded unit");

    % Now they can be safely concatenated
    T_checklist = [T_pmua T_flag];
    T_checklist = [T_pmua T_flag];
    writetable(T_checklist, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'check list (active)');
    
    % Close the notification window now that all outputs are written
    try
        if exist('h', 'var') && ishandle(h)
            fprintf('All outputs written successfully. Closing notification...\n');
            close(h);
        end
    catch ME
        fprintf('Error closing notification window: %s\n', ME.message);
    end
    fprintf('File processing completed: %s\n', fileName);
end