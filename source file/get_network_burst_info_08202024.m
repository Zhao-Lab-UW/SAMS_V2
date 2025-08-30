function get_network_burst_info_08202024(raster_raw, maxTime, samplingRate, networkParticipationThreshold, ...
    minSpikesElectrode, maxISIElectrode, minSpikesNetwork, maxISINetwork, outputFolder, sorting_results)
% GET_NETWORK_BURST_INFO_08202024 - Analyze network bursting activity across electrodes
%
% This function detects and analyzes network bursts by identifying synchronized 
% bursting activity across multiple electrodes in the recording.
%
% INPUTS:
%   rasterData - Cell array containing spike times for each electrode
%   maxTime: recordingDuration - Total recording duration in seconds
%   samplingRate - Recording sampling frequency in Hz
%   networkParticipationThreshold - Fraction of electrodes required to participate in network burst (0-1)
%   minSpikesElectrode - Minimum number of spikes per burst for electrode bursts
%   maxISIElectrode - Maximum inter-spike interval for electrode bursts (ms)
%   minSpikesNetwork - Minimum number of spikes per network burst
%   maxISINetwork - Maximum inter-spike interval for network bursts (ms)
%   outputFolder - Path to save results
%   sortingResults - Cell array containing spike sorting results
%
% OUTPUTS:
%   Excel file with network burst metrics
%   MAT file with burst information

    % Define network burst measurements
    networkMeasurements = {'number of netwwork bursts',...
        'network burst frequency (Hz)',...
        'Network burst duration -avg (s)',...
        'Network burst duration -std (s)',...
        'Number of spikes per network burst – avg',...
        'Number of spikes per network burst – std',...
        'Number of electrodes participation in burst -avg',...
        'Number of electrodes participation in burst -std',...
        'Number of spikes per network burst per channel -avg',...
        'Number of spikes per network burst per channel -std',...
        'Number of spikes per unit per network burst -avg',...
        'Number of spikes per unit per network burst -std',...
        'Network burst percentage',...
        'Network IBI coefficient of variation',...
        'Network ISI coefficient of variation',...
        'Network normalized duration IQR',...
        'Synchrony index',...
        'Number of Bursting Electrodes',...
        'Number of Active Electrodes'};

    % Initialize results table
    resultsTable = table(networkMeasurements');
    
    % Initialize storage for burst information
    burst_info_all = {};

    % Check if rasterData is empty
    if isempty(raster_raw) || size(raster_raw, 3) == 0
        fprintf('No data available for network burst analysis.\n');
        % Save empty results
        writetable(resultsTable, [outputFolder, '/spike_sorting.xlsx'], 'Sheet', 'network burst');
        return;
    end
    
    % Process each well
    for wellIndex = 1:size(raster_raw, 3)
        % Skip if no data for this well
        if isempty(raster_raw{1, 2, wellIndex})
            continue;
        end
        
        % Initialize electrode counter and data structures
        electrodeCount = 1;
        spikeTrain = zeros(size(raster_raw, 1), round(maxTime * samplingRate)); % electrodes × time
        electrodeIndices = [];
        currentElectrodeIndex = 1;
        
        % Get first electrode identifier for reference
        referencedElectrode = raster_raw{1, 2, wellIndex}(3:4);
        activeElectrodeCount = 0;
        burstingElectrodeCount = 0;
        
        % Process each electrode in this well
        for electrodeNum = 1:size(raster_raw, 1)
            % Skip empty electrodes
            if isempty(raster_raw{electrodeNum, 1, wellIndex})
                spikeTrain(end, :) = [];
            else
                % Count active electrodes
                activeElectrodeCount = activeElectrodeCount + 1;
                
                % Get spike times and convert to timepoints
                spikeTimes = round(raster_raw{electrodeNum, 1, wellIndex} * samplingRate);
                spikeTimes(spikeTimes == 0) = 1;  % Ensure valid indices
                
                % Create binary spike train
                spikeTrain(electrodeCount, spikeTimes) = 1;
                
                % Track electrode indices
                currentElectrode = raster_raw{electrodeNum, 2, wellIndex}(3:4);
                if isequal(referencedElectrode, currentElectrode)
                    electrodeIndices(electrodeCount, :) = currentElectrodeIndex;
                else
                    referencedElectrode = currentElectrode;
                    currentElectrodeIndex = currentElectrodeIndex + 1;
                    electrodeIndices(electrodeCount, :) = currentElectrodeIndex;
                end
            end
            electrodeCount = electrodeCount + 1;
        end
        
        % Only analyze if there are multiple electrodes
        if size(spikeTrain, 1) > 1
            % Detect network bursts
            networkBurstInfo = get_network_spike_participation_08072023(samplingRate, ...
                spikeTrain, electrodeIndices, networkParticipationThreshold, ...
                maxISINetwork, minSpikesNetwork);
            
            burst_info_all{wellIndex, 1} = networkBurstInfo;
            
            % Calculate network metrics
            if isempty(networkBurstInfo)
                fprintf('No network bursts detected for well %d\n', wellIndex);
                wellMetrics = zeros(19, 1);
            else
                burstingElectrodeCount = burstingElectrodeCount + 1;
                
                % Basic network burst counts and rates
                numNetworkBursts = size(networkBurstInfo, 2);
                networkBurstFrequency = numNetworkBursts / maxTime;
                
                % Duration statistics
                networkBurstDurationAvg = mean(networkBurstInfo(3, :));
                networkBurstDurationStd = std(networkBurstInfo(3, :));
                
                % Spike count statistics
                spikesPerNetworkBurstAvg = mean(networkBurstInfo(5, :));
                spikesPerNetworkBurstStd = std(networkBurstInfo(5, :));
                
                % Electrode participation statistics
                electrodesInBurstAvg = mean(networkBurstInfo(2, :));
                electrodesInBurstStd = std(networkBurstInfo(2, :));
                
                % Per-channel spike statistics
                spikesPerChannelAvg = mean(networkBurstInfo(5, :) ./ networkBurstInfo(2, :));
                spikesPerChannelStd = std(networkBurstInfo(5, :) ./ networkBurstInfo(2, :));
                
                % Per-unit spike statistics
                spikesPerUnitAvg = mean(networkBurstInfo(5, :) ./ networkBurstInfo(6, :));
                spikesPerUnitStd = std(networkBurstInfo(5, :) ./ networkBurstInfo(6, :));
                
                % Overall activity metrics
                networkBurstPercentage = sum(networkBurstInfo(5, :)) / sum(spikeTrain, 'all') * 100;
                
                % Interval variability metrics
                networkIBICoV = std(diff(networkBurstInfo(1, :))) / mean(diff(networkBurstInfo(1, :)));
                networkISICoV = std(networkBurstInfo(4, :)) / mean(networkBurstInfo(4, :));
                networkNormalizedDurationIQR = iqr(networkBurstInfo(3, :)) / median(networkBurstInfo(3, :));
                
                % Calculate synchrony index
                synchronyIndex = calculate_multivariate_synchrony(spikeTrain');
                
                % Compile all metrics
                wellMetrics = [
                    numNetworkBursts;
                    networkBurstFrequency;
                    networkBurstDurationAvg;
                    networkBurstDurationStd;
                    spikesPerNetworkBurstAvg;
                    spikesPerNetworkBurstStd;
                    electrodesInBurstAvg;
                    electrodesInBurstStd;
                    spikesPerChannelAvg;
                    spikesPerChannelStd;
                    spikesPerUnitAvg;
                    spikesPerUnitStd;
                    networkBurstPercentage;
                    networkIBICoV;
                    networkISICoV;
                    networkNormalizedDurationIQR;
                    synchronyIndex;
                    burstingElectrodeCount;
                    activeElectrodeCount
                ];
            end
        else
            % Not enough electrodes for network analysis
            wellMetrics = zeros(19, 1);
        end
        
        % Add well metrics to results table
        if ~isempty(raster_raw{1, 2, wellIndex})
            wellRow = raster_raw{1, 2, wellIndex}(1);
            wellCol = raster_raw{1, 2, wellIndex}(2);
            wellName = [char(wellRow + 'A' - 1), num2str(wellCol)];
            
            wellResultsTable = table(wellMetrics);
            wellResultsTable.Properties.VariableNames(1) = {wellName};
            resultsTable = [resultsTable, wellResultsTable];
        end
        
        % Store well identifier
        burst_info_all{wellIndex, 2} = raster_raw{1, 2, wellIndex};
    end
    
    % Save results to Excel
    writetable(resultsTable, [outputFolder, '/spike_sorting.xlsx'], 'Sheet', 'network burst');
    
    % Save detailed burst information to MAT file
    save([outputFolder, '/burst_info_all.mat'], 'burst_info_all', 'raster_raw', 'maxTime', 'sorting_results', '-v7.3');
    
    fprintf('Network burst analysis completed: %d wells processed\n', size(raster_raw, 3));
end