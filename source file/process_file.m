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

     % Start log file
    logPath = fullfile(parentFolder, 'processing_log.txt');
    logFile = fopen(logPath, 'a');
    fprintf(logFile, 'Processing file: %s\n', currentFile.name);

    % Extract file information
    fileName = currentFile.name;
    baseFileName = fileName(1:end-4); % Remove .spk extension
    
    % Create output folder for this file
    outputFolder = fullfile(parentFolder, baseFileName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    electrode_DTW_folder = fullfile(outputFolder,'DTW figures\');
    if ~exist(electrode_DTW_folder,'dir')
        mkdir(electrode_DTW_folder)
    end
    
    % Log processing
    fprintf(logFile,'Processing file %d/%d: %s\n', file_idx, total_files, fileName);
    
    % Load spike data
    filePath = fullfile(file_folder, fileName);
    try
        allData = AxisFile(filePath).SpikeData.LoadData;
        SPK_dataset = AxisFile(filePath).SpikeData;
        fr = SPK_dataset.SamplingFrequency; % Hz
        fprintf(logFile,'Sampling frequency: %d Hz\n', fr);
        
        % Get data dimensions
        [nwr, nwc, nec, ner] = size(allData);
        fprintf(logFile,'Data dimensions: [%d, %d, %d, %d]\n', nwr, nwc, nec, ner);
    catch ME
        fprintf(logFile,'Error loading file: %s\nMessage: %s\n', filePath, ME.message);
        return;
    end
    
    % Initialize analysis tables
    [T, T_electrode, T_parameters] = initialize_tables(params);
    writetable(T_parameters, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'parameter list');
    
    % Initialize data structures
    raster_raw = {};
    sorting_results = {};
    num_well = 1;
    total_num_cell = 1;
    
    % Initialize DTW results collection
    all_DTW_tables = {};
    DTW_table_count = 0;
    
    % Get recording properties
    [maxTime, maxSpikeAmp, minSpikeAmp] = get_recording_properties(allData);
    fprintf(logFile,'Recording duration: %.2f seconds\n', maxTime);
    
    % Initialize statistics
    num_total_electrodes_with_spikes = 0;
    num_inactive_units = 0;
    num_total_detected_units = 0;
    fr_inactive_units = [];
    fr_active_units = [];
    
    % Initialize PowerPoint
    try
        pptx = exportToPPTX('SAMS1.pptx');
        fprintf(logFile,'PowerPoint initialized successfully.\n');
    catch ME
        fprintf(logFile,'Error initializing PowerPoint: %s', ME.message);
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
                        fprintf(logFile,'Error opening progressbar: %s', ME.message);
                    end
                    
                    % Get current electrode data
                    Spikes = allData{i, j, m, n}(:);
                    
                    % Skip if no spikes
                    if isempty(Spikes)
                        continue;
                    end
                    
                    % Get spike times and waveforms
                    [Times, spikes] = Spikes.GetTimeVoltageVector;
                    fprintf(logFile,'Processing electrode %s: %d spikes\n', ...
                        [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)], size(spikes, 2));
                    
                    % Perform PCA
                    try
                        [~, score, ~] = pca(spikes');
                    catch ME
                        fprintf(logfile,'Error in PCA: %s\n', ME.message);
                        continue;
                    end
                    
                    % Process electrode
                    try
                        [electrode_results, electrode_stats] = process_electrode_06132025(electrode_DTW_folder, spikes, Times, i, j, m, n, ...
                            params, maxTime, maxSpikeAmp, minSpikeAmp, score);
                    catch ME
                        fprintf(logfile,'Error processing electrode: %s\n', ME.message);
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
                        fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList, ...
                        all_DTW_tables, DTW_table_count] = ...
                        update_results_06132025(electrode_results, electrode_stats, T, T_electrode, raster_raw, ...
                        sorting_results, num_electrode, total_num_cell, num_well, i, j, m, n, ...
                        num_total_electrodes_with_spikes, num_inactive_units, num_total_detected_units, ...
                        fr_inactive_units, fr_active_units, possibleMUAList, flagTooMuchList, pptx, ...
                        maxTime, maxSpikeAmp, minSpikeAmp, baseFileName, spikes, Times, score, ...
                        all_DTW_tables, DTW_table_count);
                    catch ME
                        fprintf(logfile,'Error updating results: %s\n', ME.message);
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
        fprintf(logFile,'\nFile %d/%d completed. Writing outputs...\n', file_idx, total_files);
        fprintf(logFile,'Please wait. Do not close MATLAB.\n');
    end
    
    % COMBINE AND SAVE DTW RESULTS
    if DTW_table_count > 0
        fprintf(logFile,'Combining DTW results from %d electrodes...\n', DTW_table_count);
        
        try
            % Concatenate all DTW tables
            combined_DTW_table = vertcat(all_DTW_tables{:});
            
            % Add summary statistics
            if height(combined_DTW_table) > 0
                % Calculate summary statistics
                mean_DTW_distance = mean(combined_DTW_table.DTW_Distance);
                median_DTW_distance = median(combined_DTW_table.DTW_Distance);
                std_DTW_distance = std(combined_DTW_table.DTW_Distance);
                min_DTW_distance = min(combined_DTW_table.DTW_Distance);
                max_DTW_distance = max(combined_DTW_table.DTW_Distance);
                
                % Count merge decisions
                merge_yes_count = sum(strcmp(combined_DTW_table.Merge_Decision, 'YES'));
                merge_no_count = sum(strcmp(combined_DTW_table.Merge_Decision, 'NO'));
                total_pairs = height(combined_DTW_table);
                merge_percentage = (merge_yes_count / total_pairs) * 100;
                
                % Calculate firing rate statistics
                all_firing_rates = [combined_DTW_table.Unit1_FiringRate_Hz; combined_DTW_table.Unit2_FiringRate_Hz];
                valid_firing_rates = all_firing_rates(~isnan(all_firing_rates));
                
                if ~isempty(valid_firing_rates)
                    mean_firing_rate = mean(valid_firing_rates);
                    median_firing_rate = median(valid_firing_rates);
                    std_firing_rate = std(valid_firing_rates);
                    min_firing_rate = min(valid_firing_rates);
                    max_firing_rate = max(valid_firing_rates);
                else
                    mean_firing_rate = NaN;
                    median_firing_rate = NaN;
                    std_firing_rate = NaN;
                    min_firing_rate = NaN;
                    max_firing_rate = NaN;
                end
                
                % Calculate spike count statistics
                all_spike_counts = [combined_DTW_table.Unit1_NumSpikes; combined_DTW_table.Unit2_NumSpikes];
                mean_spike_count = mean(all_spike_counts);
                median_spike_count = median(all_spike_counts);
                std_spike_count = std(all_spike_counts);
                min_spike_count = min(all_spike_counts);
                max_spike_count = max(all_spike_counts);
                
                % Create summary table
                summary_stats = {
                    'Total Unit Pairs Compared', total_pairs;
                    'Pairs Merged (YES)', merge_yes_count;
                    'Pairs Not Merged (NO)', merge_no_count;
                    'Merge Percentage (%)', merge_percentage;
                    '', '';
                    'DTW Distance Statistics', '';
                    'Mean DTW Distance', mean_DTW_distance;
                    'Median DTW Distance', median_DTW_distance;
                    'Std DTW Distance', std_DTW_distance;
                    'Min DTW Distance', min_DTW_distance;
                    'Max DTW Distance', max_DTW_distance;
                    '', '';
                    'Firing Rate Statistics (Hz)', '';
                    'Mean Firing Rate', mean_firing_rate;
                    'Median Firing Rate', median_firing_rate;
                    'Std Firing Rate', std_firing_rate;
                    'Min Firing Rate', min_firing_rate;
                    'Max Firing Rate', max_firing_rate;
                    '', '';
                    'Spike Count Statistics', '';
                    'Mean Spike Count', mean_spike_count;
                    'Median Spike Count', median_spike_count;
                    'Std Spike Count', std_spike_count;
                    'Min Spike Count', min_spike_count;
                    'Max Spike Count', max_spike_count;
                    '', '';
                    'Analysis Parameters', '';
                    'Merge Threshold Used', params.threshold_to_merge;
                    'Total Electrodes with Multiple Units', DTW_table_count;
                    'Recording Duration (sec)', maxTime
                };
                
                T_DTW_summary = table(summary_stats(:,1), summary_stats(:,2), ...
                    'VariableNames', {'Statistic', 'Value'});
                
                % Save DTW results to Excel
                fprintf(logFile,'Writing DTW results to Excel...\n');
                writetable(combined_DTW_table, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'DTW_Distances');
                writetable(T_DTW_summary, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'DTW_Summary');
                
                fprintf(logFile,'DTW analysis complete: %d unit pairs analyzed, %.1f%% merged\n', ...
                    total_pairs, merge_percentage);
                fprintf(logFile,'DTW distance range: %.4f - %.4f (mean: %.4f)\n', ...
                    min_DTW_distance, max_DTW_distance, mean_DTW_distance);
            else
                fprintf(logFile,'No DTW data to save.\n');
            end
        catch ME
            fprintf(logFile,'Error processing DTW results: %s\n', ME.message);
        end
    else
        fprintf(logFile,'No electrodes had multiple units for DTW comparison.\n');
    end
    
    % Save PowerPoint
    if ~isempty(pptx)
        try
            pptx.save([outputFolder, '\', baseFileName]);
            fprintf(logFile,'PowerPoint presentation saved successfully.\n');
        catch ME
            fprintf(logFile,'Error saving PowerPoint: %s', ME.message);
        end
    end
    
    % Save burst info
    save([outputFolder, '\burst_info_all.mat'], 'raster_raw', 'maxTime', 'sorting_results', '-v7.3');
    
    % Perform network burst analysis
    fprintf(logFile,'Performing network burst analysis...\n');
    try
        get_network_burst_info_08202024(raster_raw, maxTime, fr, params.network_participation_threshold, ...
            params.min_spikes_E, params.max_ISI_E, params.min_spikes_N, params.max_ISI_N, ...
            outputFolder, sorting_results);
    catch ME
        fprintf(logFile,'Error in network burst analysis: %s\n', ME.message);
    end
    
    % Save results to Excel
    fprintf(logFile,'Writing results to Excel...\n');
    try
        writetable(T, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'individual unit');
        writetable(T_electrode, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'electrode statistics');
    catch ME
        fprintf(logFile,'Error writing main results to Excel: %s\n', ME.message);
    end

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
    try
        T_pmua = cell2table(possibleMUAList', 'VariableNames', "Possible MultiUnit");
        T_flag = cell2table(flagTooMuchList', 'VariableNames', "over-exlcuded unit");
        T_checklist = [T_pmua T_flag];
        writetable(T_checklist, [outputFolder, '\spike_sorting.xlsx'], 'Sheet', 'check list (active)');
    catch ME
        fprintf(logFile,'Error writing check list to Excel: %s\n', ME.message);
    end
    
    % Print final summary
    fprintf(logFile,'\n=== PROCESSING SUMMARY ===\n');
    fprintf(logFile,'File: %s\n', fileName);
    fprintf(logFile,'Total electrodes with spikes: %d\n', num_total_electrodes_with_spikes);
    fprintf(logFile,'Total units detected: %d\n', num_total_detected_units);
    fprintf(logFile,'Active units: %d\n', num_total_detected_units - num_inactive_units);
    fprintf(logFile,'Inactive units: %d\n', num_inactive_units);
    if DTW_table_count > 0
        fprintf(logFile,'Electrodes with DTW analysis: %d\n', DTW_table_count);
    end
    fprintf(logFile,'Possible MUA electrodes: %d\n', length(possibleMUAList) - (strcmp(possibleMUAList{end}, 'no check list') || strcmp(possibleMUAList{end}, '')));
    fprintf(logFile,'Over-excluded electrodes: %d\n', length(flagTooMuchList) - (strcmp(flagTooMuchList{end}, 'no check list') || strcmp(flagTooMuchList{end}, '')));
    fprintf(logFile,'Recording duration: %.2f seconds\n', maxTime);
    fprintf(logFile,'Output folder: %s\n', outputFolder);
    fprintf(logFile,'=========================\n\n');
    
    % Close the notification window now that all outputs are written
    try
        if exist('h', 'var') && ishandle(h)
            fprintf(logFile,'All outputs written successfully. Closing notification...\n');
            close(h);
        end
    catch ME
        fprintf(logFile,'Error closing notification window: %s\n', ME.message);
    end
    
    fprintf(logFile,'File processing completed: %s\n', fileName);
    % Close log file
    fclose(logFile);
end