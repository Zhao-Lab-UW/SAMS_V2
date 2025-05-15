function burstInfo = get_burst(spikeTimes, maxInterSpikeInterval, minSpikesPerBurst)
% GET_BURST - Detect bursts in spike trains based on inter-spike intervals
%
% This function detects bursts in spike times based on maximum inter-spike interval
% and minimum number of spikes per burst criteria.
%
% INPUTS:
%   spikeTimes - Vector of spike times
%   maxInterSpikeInterval - Maximum time between consecutive spikes to be considered in the same burst
%   minSpikesPerBurst - Minimum number of spikes required to form a burst
%
% OUTPUTS:
%   burstInfo - Matrix containing burst information:
%       Row 1: Starting spike index of each burst
%       Row 2: Number of spikes in each burst
%       Row 3: Duration of each burst
%       Row 4: Mean ISI within each burst

    % Initialize burst detection
    isBurstSpike = false(size(spikeTimes));
    isBurstSpike(1) = true;  % First spike starts a potential burst
    
    % Identify spikes that could be part of a burst based on ISI
    for spikeIdx = 2:length(spikeTimes)
        if (spikeTimes(spikeIdx) - spikeTimes(spikeIdx-1)) < maxInterSpikeInterval
            isBurstSpike(spikeIdx) = true;
        end
    end
    
    % Initialize burst statistics
    burstCount = 0;
    burstStartIndices = [];
    spikesPerBurst = [];
    burstDurations = [];
    meanISIWithinBurst = [];
    
    % Find bursts by iterating through spike times
    spikeIdx = 1;
    while spikeIdx <= (length(spikeTimes) - (minSpikesPerBurst - 1))
        spikeCountInBurst = 0;
        
        % Check if current spike could start a burst
        if isBurstSpike(spikeIdx)
            % Count this spike
            spikeCountInBurst = 1;
            
            % Count consecutive burst spikes
            for nextSpikeIdx = (spikeIdx + 1):length(spikeTimes)
                if isBurstSpike(nextSpikeIdx)
                    spikeCountInBurst = spikeCountInBurst + 1;
                else
                    break;
                end
            end
        end
        
        % Check if we found a valid burst
        if spikeCountInBurst >= minSpikesPerBurst
            burstCount = burstCount + 1;
            
            % Store burst start index
            burstStartIndices(burstCount) = spikeIdx;
            
            % Store number of spikes in this burst
            spikesPerBurst(burstCount) = spikeCountInBurst;
            
            % Calculate burst duration
            burstDurations(burstCount) = spikeTimes(spikeIdx + spikeCountInBurst - 1) - spikeTimes(spikeIdx);
            
            % Calculate mean ISI within burst
            burstSpikeIndices = spikeIdx:(spikeIdx + spikeCountInBurst - 1);
            burstISIs = diff(spikeTimes(burstSpikeIndices));
            meanISIWithinBurst(burstCount) = mean(burstISIs);
            
            % Move to the spike after this burst
            spikeIdx = spikeIdx + spikeCountInBurst;
        else
            % Move to next spike
            spikeIdx = spikeIdx + 1;
        end
    end
    
    % Combine burst information
    burstInfo = [burstStartIndices; spikesPerBurst; burstDurations; meanISIWithinBurst];
end