function [meanDistance, mergeGroups, mergeDistances, nonMergeDistances] = template_comparison(waveforms, clusterIndices, mergeThreshold,overlap_threshold)
% TEMPLATE_COMPARISON - Compare spike templates to identify similar units for merging
%
% This function calculates distances between mean templates of spike clusters
% and identifies which clusters should be merged based on similarity.
%
% INPUTS:
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%   clusterIndices - Cell array containing spike indices for each cluster
%   mergeThreshold - Distance threshold below which clusters are merged
%
% OUTPUTS:
%   meanDistance - Mean distance between all template pairs
%   mergeGroups - Cell array of cluster indices to be merged
%   mergeDistances - Distances between templates that should be merged
%   nonMergeDistances - Distances between templates that should not be merged

    % Calculate mean template for each cluster
    templateCount = 0;
    meanTemplates = [];
    validClusterIndices = [];


    clusterIndices = clusterIndices(~cellfun('isempty',clusterIndices));
    for clusterIdx = 1:length(clusterIndices)
        templateCount = templateCount + 1;
        validClusterIndices(templateCount) = clusterIdx;
        meanTemplates(:, templateCount) = mean(waveforms(:, clusterIndices{clusterIdx}), 2);
    end
    
    % Normalize templates for better comparison
    templateMin = min(meanTemplates(:));
    templateMax = max(meanTemplates(:));
    normalizedTemplates = (meanTemplates - templateMin) / (templateMax - templateMin);
    
    % Calculate distances between templates
    pairwiseDistances = [];
    mergeCandidate = [];
    mergeCount = 0;
    nonMergeCount = 0;
    mergeDistances = [];
    nonMergeDistances = [];
    
    % Compare each template pair
    for template1_idx = 1:templateCount-1
        for template2_idx = template1_idx+1:templateCount
            % Calculate DTW distance
            currentDistance = dtw(normalizedTemplates(:, template1_idx), normalizedTemplates(:, template2_idx));    
            pairwiseDistances(end+1) = currentDistance;

            % Check if templates have similar shape features
            template1 = normalizedTemplates(:, template1_idx);
            template2 = normalizedTemplates(:, template2_idx);
            
            % Find peaks and troughs
            peak1 = max(template1);
            trough1 = min(template1);
            peak2 = max(template2);
            trough2 = min(template2);
            
            % Calculate height
            height1 = peak1 - trough1;
            height2 = peak2 - trough2;
            
            % Define threshold for peak/trough matching
            shapeSimilarityThreshold = 1-overlap_threshold; %using params.overlap_threshold=0.7, 
            
            % Check if peaks match
            peakMatch = (peak1 < (peak2 + height1 * shapeSimilarityThreshold)) && ...
                       (peak1 > (peak2 - height1 * shapeSimilarityThreshold));
            
            % Check if troughs match
            troughMatch = (trough1 < (trough2 + height1 * shapeSimilarityThreshold)) && ...
                         (trough1 > (trough2 - height1 * shapeSimilarityThreshold));
            
            % Determine if templates should be merged
            shouldMerge = (currentDistance < mergeThreshold) && peakMatch && troughMatch;
            if shouldMerge
                mergeCount = mergeCount + 1;
                mergeCandidate(mergeCount, :) = [template1_idx, template2_idx];
                mergeDistances(mergeCount) = currentDistance;
            else
                nonMergeCount = nonMergeCount + 1;
                nonMergeDistances(nonMergeCount) = currentDistance;
            end
        end
    end

    % Calculate mean distance
    if ~isempty(pairwiseDistances)
        meanDistance = mean(pairwiseDistances);
    else
        meanDistance = 999; % Default high value if no comparisons
    end
    
    % Group merge candidates into connected components
    mergeGroups = {};
    
    % Process each merge candidate pair
    for candidateIdx = 1:size(mergeCandidate, 1)
        currentPair = validClusterIndices(mergeCandidate(candidateIdx, :));
        
        if candidateIdx == 1
            % First pair starts a new group
            mergeGroups{1} = currentPair;
        else
            % Check if pair connects to existing group
            foundMatch = false;
            
            for groupIdx = 1:length(mergeGroups)
                if any(ismember(currentPair, mergeGroups{groupIdx}))
                    % Add to existing group
                    mergeGroups{groupIdx} = unique([mergeGroups{groupIdx}, currentPair]);
                    foundMatch = true;
                    break;
                end
            end
            
            if ~foundMatch
                % Create new group
                mergeGroups{end+1} = currentPair;
            end
        end
    end
end