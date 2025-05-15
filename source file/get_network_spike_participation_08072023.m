function [networkBurstInfo] = get_network_spike_participation_08072023(samplingRate, spikeMatrix, electrodeIndices, networkThreshold, maxNetworkISI, minNetworkSpikes)
% GET_NETWORK_SPIKE_PARTICIPATION_08072023 - Detect network bursts across multiple electrodes
%
% This function identifies synchronized bursting across multiple electrodes
% based on network participation threshold and burst criteria.
%
% INPUTS:
%   samplingRate - Recording sampling rate in Hz
%   spikeMatrix - Binary spike matrix (electrodes × timepoints)
%   electrodeIndices - Indices of electrodes in recording
%   networkThreshold - Minimum fraction of electrodes that must participate
%   maxNetworkISI - Maximum inter-spike interval for network burst (ms)
%   minNetworkSpikes - Minimum spikes required for a network burst
%
% OUTPUTS:
%   networkBurstInfo - Matrix containing network burst information:
%       Row 1: Starting timepoint of each burst
%       Row 2: Number of electrodes participating in each burst
%       Row 3: Duration of each burst (seconds)
%       Row 4: Mean ISI within each burst (seconds)
%       Row 5: Number of spikes per network burst
%       Row 6: Number of firing cells in each burst

    % Convert maxNetworkISI from ms to timepoints
    maxNetworkISI_timepoints = maxNetworkISI * samplingRate / 1000;
    
    % Find timepoints with any firing
    totalFiringByTimepoint = sum(spikeMatrix);
    firingTimepoints = find(totalFiringByTimepoint > 0);
    
    % Initialize burst detection
    isBurstTimepoint = false(size(firingTimepoints));
    isBurstTimepoint(1) = true;  % First timepoint starts a potential burst
    
    % Find timepoints that could be part of a burst based on ISI
    for timepointIdx = 2:length(firingTimepoints)
        if (firingTimepoints(timepointIdx) - firingTimepoints(timepointIdx-1)) < maxNetworkISI_timepoints
            isBurstTimepoint(timepointIdx) = true;
        end
    end
    
    % Initialize network burst statistics
    burstCount = 0;
    burstStartTimepoints = [];
    electrodesPerBurst = [];
    burstDurations = [];
    meanISIWithinBurst = [];
    spikesPerBurst = [];
    cellsPerBurst = [];
    
    % Find network bursts
    timepointIdx = 1;
    while timepointIdx <= length(firingTimepoints)
        if isBurstTimepoint(timepointIdx)
            % Start tracking a potential burst
            burstStartIdx = timepointIdx;
            firingElectrodeIndices = [];
            firingCellIndices = [];
            
            % Find indices of firing electrodes at the start
            currentFiringElectrodes = find(spikeMatrix(:, firingTimepoints(timepointIdx)) == 1);
            firingElectrodeIndices = [firingElectrodeIndices; electrodeIndices(currentFiringElectrodes)];
            firingCellIndices = [firingCellIndices; find(spikeMatrix(:, firingTimepoints(timepointIdx)) == 1)];
            
            % Count consecutive burst timepoints
            consecutiveCount = 1;
            for nextTimepointIdx = (timepointIdx + 1):length(firingTimepoints)
                if isBurstTimepoint(nextTimepointIdx)
                    consecutiveCount = consecutiveCount + 1;
                    
                    % Add new firing electrodes to list
                    currentFiringElectrodes = find(spikeMatrix(:, firingTimepoints(nextTimepointIdx)) == 1);
                    firingElectrodeIndices = [firingElectrodeIndices; electrodeIndices(currentFiringElectrodes)];
                    firingCellIndices = [firingCellIndices; find(spikeMatrix(:, firingTimepoints(nextTimepointIdx)) == 1)];
                else
                    break;
                end
            end
            
            % Get unique electrode indices
            uniqueElectrodeIndices = unique(firingElectrodeIndices);
            uniqueCellIndices = unique(firingCellIndices);
            
            % Check if enough electrodes participated and enough spikes occurred
            totalElectrodes = length(unique(electrodeIndices));
            totalSpikes = sum(spikeMatrix(:, firingTimepoints(burstStartIdx:(burstStartIdx+consecutiveCount-1))), 'all');
            
            if (length(uniqueElectrodeIndices) > networkThreshold * totalElectrodes) && ...
               (totalSpikes >= minNetworkSpikes)
                
                % Valid network burst detected
                burstCount = burstCount + 1;
                
                % Store burst information
                burstStartTimepoints(burstCount) = firingTimepoints(burstStartIdx);
                electrodesPerBurst(burstCount) = length(uniqueElectrodeIndices);
                burstDurations(burstCount) = firingTimepoints(burstStartIdx + consecutiveCount - 1) - firingTimepoints(burstStartIdx);
                
                % Calculate mean ISI
                burstTimepointIndices = burstStartIdx:(burstStartIdx + consecutiveCount - 1);
                burstTimepoints = firingTimepoints(burstTimepointIndices);
                meanISIWithinBurst(burstCount) = mean(diff(burstTimepoints));
                
                % Store spike and cell counts
                spikesPerBurst(burstCount) = totalSpikes;
                cellsPerBurst(burstCount) = length(uniqueCellIndices);
                
                % Skip to after this burst
                timepointIdx = burstStartIdx + consecutiveCount;
            else
                % Not a valid network burst, move to next timepoint
                timepointIdx = timepointIdx + 1;
            end
        else
            % Move to next timepoint
            timepointIdx = timepointIdx + 1;
        end
    end
    
    % Convert durations and ISIs from timepoints to seconds
    burstDurations = burstDurations / samplingRate;
    meanISIWithinBurst = meanISIWithinBurst / samplingRate;
    
    % Combine network burst information
    networkBurstInfo = [burstStartTimepoints; electrodesPerBurst; burstDurations; ...
                         meanISIWithinBurst; spikesPerBurst; cellsPerBurst];
    
    % Replace any NaN values with zeros
    networkBurstInfo(isnan(networkBurstInfo)) = 0;
end