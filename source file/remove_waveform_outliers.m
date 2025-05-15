
function cleanWaveforms = remove_waveform_outliers(waveforms)
% REMOVE_WAVEFORM_OUTLIERS - Remove outlier waveforms based on amplitude
%
% INPUTS:
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%
% OUTPUTS:
%   cleanWaveforms - Matrix of waveforms with outliers removed

    cleanWaveforms = waveforms;
    
    % Check each time point for outliers
    for timeIdx = 1:size(waveforms, 1)
        [~, outlierMask] = rmoutliers(waveforms(timeIdx, :), 'mean');
        outlierIndices = find(outlierMask == 1);
        
        % Remove outliers from all time points
        if ~isempty(outlierIndices)
            cleanWaveforms(:, outlierIndices) = [];
        end
    end
end