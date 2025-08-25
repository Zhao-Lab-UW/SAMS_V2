function [maxTime, maxSpikeAmp, minSpikeAmp] = get_recording_properties(allData)
    % Get recording duration, max and min spike amplitudes
    maxTime = 0;
    maxSpikeAmps = [];
    minSpikeAmps = [];
    [nwr, nwc, nec, ner] = size(allData);
    
    for i = 1:nwr
        for j = 1:nwc
            for m = 1:nec
                for n = 1:ner
                    spikes = allData{i,j,m,n}(:);
                    if ~isempty(spikes)
                        [Times, spikes] = spikes.GetTimeVoltageVector;
                        tmp = max(Times(:));
                        maxSpikeAmps = [maxSpikeAmps max(spikes)];
                        minSpikeAmps = [minSpikeAmps min(spikes)];
                        if maxTime < tmp
                            maxTime = tmp;
                        end
                    end
                end
            end
        end
    end
    
    % Get 95% maxSpikeAmps for normalization
    maxSpikeAmp = sort(maxSpikeAmps);
    maxSpikeAmp = maxSpikeAmp(round(length(maxSpikeAmps)*0.95));
    
    minSpikeAmp = sort(minSpikeAmps);
    minSpikeAmp = minSpikeAmp(round(length(minSpikeAmps)*0.05)+1);
    
    fprintf('Recording duration: %.2f seconds\n', maxTime);
    fprintf('95th percentile spike amplitude: %.4f\n', maxSpikeAmp);
end