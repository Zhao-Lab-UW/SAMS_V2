
function cleanWaveforms = remove_waveform_outliers(waveforms)
% REMOVE_WAVEFORM_OUTLIERS - Remove outlier waveforms based on amplitude
%
% INPUTS:
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%
% OUTPUTS:
%   cleanWaveforms - Matrix of waveforms with outliers removed

cleanWaveforms = waveforms;

% Safe approach: accumulate outliers, remove once
allOutliers = false(1, size(waveforms, 2));

for timeIdx = 1:size(waveforms, 1)
    [~, outlierMask] = rmoutliers(waveforms(timeIdx, :), 'mean'); % Use original matrix
    allOutliers = allOutliers | outlierMask; % Accumulate outliers
end

cleanWaveforms = waveforms(:, ~allOutliers); % Remove all at once
end