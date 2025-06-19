function [T_unit, T_electrode] = get_individual_unit_analysis(sorting_results, allData, maxTime, ...
    cutoffFrequency, maxISI, minSpikesPerBurst, outputFolder)
% GET_INDIVIDUAL_UNIT_ANALYSIS - Extract and analyze individual unit statistics
%
% This function loops through sorting_results to analyze individual units and electrodes,
% extracting spike times, calculating metrics, and generating output tables.
%
% INPUTS:
%   sorting_results - Cell array containing sorted unit indices
%   allData - Cell array containing spike data from all electrodes
%   maxTime - Maximum recording duration in seconds
%   cutoffFrequency - Cutoff frequency for filtering (Hz) - units with firing rates 
%                     below this value will be excluded from the analysis
%   maxISI - Maximum inter-spike interval for burst detection (ms)
%   minSpikesPerBurst - Minimum spikes per burst
%   outputFolder - Folder to save results
%
% OUTPUTS:
%   T_unit - Table with individual unit statistics
%   T_electrode - Table with electrode-level statistics

    % Define the measurements to track for individual units
    unitMeasurements = {'Number of Spikes';
        'Mean Firing Rate (Hz)';
        'ISI in general';
        'Number of Bursts';
        'Burst Frequency (Hz)';
        'Burst Duration -Avg(sec)';
        'Burst Duration -Std(sec)';
        'Burst Duration -Median(sec)';
        'Burst Percentage (%)';
        'Normalized Burst Duration IQR';
        'Number of Spikes per Burst -Avg';
        'Number of Spikes per Burst -Std';
        'Number of Spikes per Burst -Median';
        'Mean ISI within Burst -Avg';
        'Mean ISI within Burst -Std';
        'Mean ISI within Burst -Median';
        'ISI Coefficient of Variation';
        'Inter-Burst Interval - Avg';
        'Inter-Burst Interval - Std';
        'Inter-Burst Interval - Median';
        'Inter-Burst Interval Coefficient of Variation';
        'Amplitude -Avg(mV)';
        'Amplitude -Std(mV)';
        'Amplitude - Avg Peak to Trough(mV)';
        'Amplitude - Std Peak to Trough(mV)';
        'max ISI (ms)';
        'number of removed short ISIs'};
    
    % Define electrode-level measurements
    electrodeMeasurements = {
        'Total number of spikes';
        'Number of excluded spikes';
        'Unsorted firing rates (Hz)';
        'Number of Units Detected';
        'Number of Active Units (>= cutoff)';
        'Number of Inactive Units (< cutoff)'
    };

    % Define colors for units - extended to support 10 units
    colorList = {
        [1, 0, 0],       % red
        [0, 0.7, 0],     % green
        [0, 0, 1],       % blue
        [0, 0.8, 0.8],   % cyan
        [1, 0, 1],       % magenta
        [0.8, 0.8, 0],   % yellow
        [0, 0, 0],       % black
        [1, 0.5, 0],     % orange
        [0.5, 0, 0.5],   % purple
        [0, 0.5, 0]      % dark green
    };
    
    % Initialize output matrices instead of tables first
    unitData = zeros(length(unitMeasurements), 0);
    electrodeData = zeros(length(electrodeMeasurements), 0);
    
    % Arrays to store column names
    unitNames = {};
    electrodeNames = {};
    
    % Get dimensions of sorting_results array
    [nwr, nwc, nec, ner, ~] = size(sorting_results);
    
    % Convert maxISI from ms to seconds for internal calculations
    maxISI_seconds = maxISI / 1000;
    
    % Loop through all electrodes
    for i = 1:nwr
        for j = 1:nwc
            for m = 1:nec
                for n = 1:ner
                    % Skip if no sorting results for this electrode
                    if isempty(sorting_results{i,j,m,n,1})
                        continue;
                    end
                    
                    % Get electrode name
                    electrodeName = [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)];
                    
                    % Get spike data for this electrode
                    Spikes = allData{i,j,m,n}(:);
                    if isempty(Spikes)
                        continue;
                    end
                    
                    % Get spike times and waveforms
                    [Times, spikes] = Spikes.GetTimeVoltageVector;
                    
                    % Get sorted units and unsorted waveforms
                    sorted_units = sorting_results{i,j,m,n,1};
                    unsorted_wfs = sorting_results{i,j,m,n,2};
                    
                    % Calculate electrode-level statistics
                    total_spikes = size(spikes, 2);
                    excluded_spikes = length(unsorted_wfs);
                    unsorted_fr = excluded_spikes / maxTime;
                    num_units = length(sorted_units);
                    
                    % Count active vs inactive units based on firing rate
                    active_units = 0;
                    inactive_units = 0;
                    
                    for unit_idx = 1:length(sorted_units)
                        if isempty(sorted_units{unit_idx})
                            continue;
                        end
                        
                        % Calculate firing rate for this unit
                        unit_fr = length(sorted_units{unit_idx}) / maxTime;
                        
                        if unit_fr >= cutoffFrequency
                            active_units = active_units + 1;
                        else
                            inactive_units = inactive_units + 1;
                        end
                    end
                    
                    % Add electrode statistics to data matrix
                    electrodeNames{end+1} = electrodeName;
                    electrodeData(:, end+1) = [
                        total_spikes;
                        excluded_spikes;
                        unsorted_fr;
                        num_units;
                        active_units;
                        inactive_units
                    ];
                    
                    % Process each unit for this electrode
                    for unit_idx = 1:length(sorted_units)
                        if isempty(sorted_units{unit_idx})
                            continue;
                        end
                        
                        % Get unit name
                        unitName = [electrodeName, '_', num2str(unit_idx)];
                        
                        % Get spike indices for this unit
                        unit_spike_indices = sorted_units{unit_idx};
                        
                        % Get spike times and waveforms for this unit
                        unit_times = Times(unit_spike_indices);
                        unit_waveforms = spikes(:, unit_spike_indices);
                        
                        % Number of spikes
                        num_spikes = length(unit_times);
                        
                        % Mean firing rate
                        mean_fr = num_spikes / maxTime;
                        
                        % Skip this unit if firing rate is below cutoff frequency
                        if mean_fr < cutoffFrequency
                            fprintf('Unit %s firing rate (%.2f Hz) below cutoff (%.2f Hz) - excluding from analysis\n', ...
                                unitName, mean_fr, cutoffFrequency);
                            continue;
                        end
                        
                        % Calculate ISIs
                        if length(unit_times) > 1
                            isis = diff(unit_times);
                            mean_isi = mean(isis);
                            isi_cv = std(isis) / mean_isi;  % Coefficient of variation
                        else
                            isis = [];
                            mean_isi = NaN;
                            isi_cv = NaN;
                        end
                        
                        % Detect bursts using the get_burst function
                        if length(unit_times) >= minSpikesPerBurst
                            burstInfo = get_burst(unit_times, maxISI_seconds, minSpikesPerBurst);
                        else
                            burstInfo = [];
                        end
                        
                        % Calculate burst statistics
                        if ~isempty(burstInfo) && size(burstInfo, 2) > 0
                            num_bursts = size(burstInfo, 2);
                            burst_freq = num_bursts / maxTime;
                            
                            % Burst duration stats
                            burst_durations = burstInfo(3, :);
                            burst_duration_avg = mean(burst_durations);
                            burst_duration_std = std(burst_durations);
                            burst_duration_median = median(burst_durations);
                            burst_duration_iqr = iqr(burst_durations);
                            norm_burst_duration_iqr = burst_duration_iqr / burst_duration_median;
                            
                            % Spikes per burst stats
                            spikes_per_burst = burstInfo(2, :);
                            spikes_per_burst_avg = mean(spikes_per_burst);
                            spikes_per_burst_std = std(spikes_per_burst);
                            spikes_per_burst_median = median(spikes_per_burst);
                            
                            % ISI within burst stats
                            isi_within_burst = burstInfo(4, :);
                            isi_within_burst_avg = mean(isi_within_burst);
                            isi_within_burst_std = std(isi_within_burst);
                            isi_within_burst_median = median(isi_within_burst);
                            
                            % Inter-burst interval stats
                            if num_bursts > 1
                                burst_start_times = unit_times(burstInfo(1, :));
                                ibis = diff(burst_start_times);
                                ibi_avg = mean(ibis);
                                ibi_std = std(ibis);
                                ibi_median = median(ibis);
                                ibi_cv = ibi_std / ibi_avg;
                            else
                                ibi_avg = NaN;
                                ibi_std = NaN;
                                ibi_median = NaN;
                                ibi_cv = NaN;
                            end
                            
                            % Calculate burst percentage (% of spikes in bursts)
                            total_spikes_in_bursts = sum(spikes_per_burst);
                            burst_percentage = (total_spikes_in_bursts / num_spikes) * 100;
                        else
                            % No bursts detected
                            num_bursts = 0;
                            burst_freq = 0;
                            burst_duration_avg = NaN;
                            burst_duration_std = NaN;
                            burst_duration_median = NaN;
                            norm_burst_duration_iqr = NaN;
                            spikes_per_burst_avg = NaN;
                            spikes_per_burst_std = NaN;
                            spikes_per_burst_median = NaN;
                            isi_within_burst_avg = NaN;
                            isi_within_burst_std = NaN;
                            isi_within_burst_median = NaN;
                            ibi_avg = NaN;
                            ibi_std = NaN;
                            ibi_median = NaN;
                            ibi_cv = NaN;
                            burst_percentage = 0;
                        end
                        
                        % Calculate waveform amplitude statistics
                        if ~isempty(unit_waveforms)
                            % Average amplitude
                            amp_avg = mean(max(unit_waveforms) - min(unit_waveforms));
                            amp_std = std(max(unit_waveforms) - min(unit_waveforms));
                            
                            % Peak to trough
                            peak_to_trough = max(unit_waveforms) - min(unit_waveforms);
                            peak_to_trough_avg = mean(peak_to_trough);
                            peak_to_trough_std = std(peak_to_trough);
                        else
                            amp_avg = NaN;
                            amp_std = NaN;
                            peak_to_trough_avg = NaN;
                            peak_to_trough_std = NaN;
                        end
                        
                        % Count removed short ISIs (ISIs < refractory period)
                        % Assuming refractory period is 1.5 ms
                        if ~isempty(isis)
                            removed_short_isis = sum(isis < 0.0015);
                        else
                            removed_short_isis = 0;
                        end
                        
                        % Compile unit statistics and add to data matrix
                        unitNames{end+1} = unitName;
                        unitData(:, end+1) = [
                            num_spikes;
                            mean_fr;
                            mean_isi;
                            num_bursts;
                            burst_freq;
                            burst_duration_avg;
                            burst_duration_std;
                            burst_duration_median;
                            burst_percentage;
                            norm_burst_duration_iqr;
                            spikes_per_burst_avg;
                            spikes_per_burst_std;
                            spikes_per_burst_median;
                            isi_within_burst_avg;
                            isi_within_burst_std;
                            isi_within_burst_median;
                            isi_cv;
                            ibi_avg;
                            ibi_std;
                            ibi_median;
                            ibi_cv;
                            amp_avg;
                            amp_std;
                            peak_to_trough_avg;
                            peak_to_trough_std;
                            maxISI;
                            removed_short_isis;
                        ];
                    end
                end
            end
        end
    end
    
    % Convert data matrices to tables
    if isempty(unitData)
        % Create empty table with correct variable names if no units pass the filter
        T_unit = array2table(zeros(length(unitMeasurements), 0), 'RowNames', unitMeasurements);
    else
        % Create table with unit data
        T_unit = array2table(unitData, 'RowNames', unitMeasurements, 'VariableNames', unitNames);
    end
    
    if isempty(electrodeData)
        % Create empty table with correct variable names if no electrodes
        T_electrode = array2table(zeros(length(electrodeMeasurements), 0), 'RowNames', electrodeMeasurements);
    else
        % Create table with electrode data
        T_electrode = array2table(electrodeData, 'RowNames', electrodeMeasurements, 'VariableNames', electrodeNames);
    end
    
    % Save intermediate results
    save(fullfile(outputFolder, 'unit_analysis.mat'), 'T_unit', 'T_electrode', '-v7.3');
end