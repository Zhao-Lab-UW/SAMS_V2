function [mergedOutlierList, unsortedWaveformClusters] = template_comparison_outlier_02212024(waveforms, templateClusters, outlierIndices, mergeThreshold, clusterFunction, kValues, pcaScores, stdCutoff)
% TEMPLATE_COMPARISON_OUTLIER_02212024 - Compare outlier waveforms with existing templates
%
% This function compares outlier waveforms with existing templates to determine if
% they can be merged with existing clusters. It also attempts to cluster remaining
% outliers into new potential units.
%
% INPUTS:
%   waveforms - Spike waveforms matrix (timepoints × spikes)
%   templateClusters - Cell array of indices for each existing cluster
%   outlierIndices - Indices of outlier waveforms
%   mergeThreshold - Distance threshold below which waveforms are merged
%   clusterFunction - Function handle for clustering algorithm
%   kValues - Range of k values to try for clustering
%   pcaScores - PCA scores for all waveforms
%   stdCutoff - Standard deviation cutoff for outlier removal
%
% OUTPUTS:
%   mergedOutlierList - Cell array of outlier indices to merge with each cluster
%   unsortedWaveformClusters - Cell array of new clusters from remaining outliers

    % Calculate mean waveform templates for each cluster
    meanTemplates = [];
    clusterCount = 0;
    
    for templateIndex = 1:length(templateClusters)
        if ~isempty(templateClusters{templateIndex})
            clusterCount = clusterCount + 1;
            meanTemplates(:, clusterCount) = mean(waveforms(:, templateClusters{templateIndex}), 2);
        end
    end

    % Normalize templates for better comparison
    templateMin = min(meanTemplates(:));
    templateMax = max(meanTemplates(:));
    normalizedTemplates = (meanTemplates - templateMin) / (templateMax - templateMin);
    
    % Normalize waveforms for comparison
    normalizedWaveforms = (waveforms - min(waveforms(:))) / (max(waveforms(:)) - min(waveforms(:)));
    normalizedWaveforms = normalizedWaveforms(1:size(normalizedWaveforms, 1), :);

    % Initialize variables for tracking distances and matches
    mergeList = [];
    mergeCount = 0;
    distancesForOutliers = [];
    
    % Compare each outlier waveform with each template
    for outlierIndex = 1:length(outlierIndices)
        distances = [];
        matchFlags = [];
        outlierWaveform = normalizedWaveforms(:, outlierIndices(outlierIndex));
        
        % Compare with each template
        for templateIndex = 1:clusterCount
            % Calculate DTW distance
            distance = dtw(outlierWaveform, normalizedTemplates(:, templateIndex));
            distances(templateIndex) = distance;
            
            % Check feature similarity
            outlierPeak = max(outlierWaveform);
            outlierTrough = min(outlierWaveform);
            templatePeak = max(normalizedTemplates(:, templateIndex));
            templateTrough = min(normalizedTemplates(:, templateIndex));
            
            outlierHeight = outlierPeak - outlierTrough;
            templateHeight = templatePeak - templateTrough;
            
            shapeThreshold = 0.3;
            peakMatch = (outlierPeak < (templatePeak + outlierHeight * shapeThreshold)) && ...
                        (outlierPeak > (templatePeak - outlierHeight * shapeThreshold));
            troughMatch = (outlierTrough < (templateTrough + outlierHeight * shapeThreshold)) && ...
                          (outlierTrough > (templateTrough - outlierHeight * shapeThreshold));
            
            % Record if template and outlier match in shape features
            if peakMatch && troughMatch
                matchFlags(templateIndex) = 1;
            else
                matchFlags(templateIndex) = 0;
            end
        end
        
        % Find best matching template
        [minDistance, bestTemplateIndex] = min(distances);
        
        % Check if this outlier should be merged with a template
        if matchFlags(bestTemplateIndex) == 1
            mergeCount = mergeCount + 1;
            mergeList(mergeCount, :) = [outlierIndex, bestTemplateIndex];
        end
        
        % Store distances for debugging/analysis
        distancesForOutliers(outlierIndex, :) = distances;
    end
    
    % Create list of outliers to merge with each template
    mergedOutlierList = cell(1, length(templateClusters));
    mergedOutliers = [];
    
    if ~isnan(mergeList)
        for templateIndex = 1:length(templateClusters)
            % Find outliers that match this template
            matchIndices = find(mergeList(:, 2) == templateIndex);
            
            % Add their indices to the merge list
            mergedOutlierList{templateIndex} = outlierIndices(mergeList(matchIndices, 1));
            
            % Track which outliers have been merged
            mergedOutliers = [mergedOutliers; outlierIndices(mergeList(matchIndices, 1))];
        end
    end
    
    % Process remaining outliers that weren't merged
    unsortedOutliers = setdiff(outlierIndices, mergedOutliers);
    unsortedWaveformClusters = {};
    
    % Attempt to cluster remaining outliers
    if ~isempty(unsortedOutliers)
        % Extract features for clustering (peaks and valleys)
        peakValues = max(waveforms(:, unsortedOutliers));
        troughValues = min(waveforms(:, unsortedOutliers));
        features = [peakValues; troughValues]';
        
        % Use PCA scores if available, otherwise use peak/trough features
        if ~isempty(pcaScores)
            try
                % Evaluate optimal clustering
                evalResult = evalclusters(pcaScores(unsortedOutliers, 1:2), clusterFunction, "DaviesBouldin", 'klist', kValues);
                
                if ~isnan(evalResult.OptimalK)
                    % Perform clustering
                    [clusterIndices, ~] = spectralcluster(pcaScores(unsortedOutliers, 1:2), evalResult.OptimalK);
                    
                    % Process each new cluster
                    for newClusterIndex = 1:max(clusterIndices)
                        % Get indices for this cluster
                        clusterMemberIndices = unsortedOutliers(clusterIndices == newClusterIndex);
                        
                        % Get waveforms for this cluster
                        clusterWaveforms = waveforms(:, clusterMemberIndices);
                        
                        % Remove outliers within this cluster
                        for timeIndex = 1:size(clusterWaveforms, 1)
                            [~, outlierFlags] = rmoutliers(clusterWaveforms(timeIndex, :), 'mean', 'ThresholdFactor', stdCutoff);
                            clusterMemberIndices(outlierFlags) = [];
                        end
                        
                        % Store cleaned cluster
                        unsortedWaveformClusters{newClusterIndex} = clusterMemberIndices;
                    end
                end
            catch
                % Fallback to simpler clustering if PCA-based clustering fails
                unsortedWaveformClusters = {};
            end
        end
    end
end