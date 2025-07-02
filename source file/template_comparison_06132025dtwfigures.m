function [meanDistance, mergeGroups, mergeDistances, nonMergeDistances, T_electrode_DTW] = template_comparison_06132025dtwfigures(electrode_DTW_folder,i,j,m,n,Times,waveforms, clusterIndices, mergeThreshold)
% TEMPLATE_COMPARISON_06122025DTWFIGURES - Compare spike templates to identify similar units for merging
%
% This function calculates distances between mean templates of spike clusters
% and identifies which clusters should be merged based on similarity.
%
% INPUTS:
%   electrode_DTW_folder - Folder to save DTW figures
%   i,j,m,n - Electrode indices
%   Times - Spike timing vector
%   waveforms - Matrix of spike waveforms (timepoints × spikes)
%   clusterIndices - Cell array containing spike indices for each cluster
%   mergeThreshold - Distance threshold below which clusters are merged
%
% OUTPUTS:
%   meanDistance - Mean distance between all template pairs
%   mergeGroups - Cell array of cluster indices to be merged
%   mergeDistances - Distances between templates that should be merged
%   nonMergeDistances - Distances between templates that should not be merged
%   T_electrode_DTW - Table with DTW distances and unit statistics

    % Calculate mean template for each cluster
    templateCount = 0;
    meanTemplates = [];
    validClusterIndices = [];
    currElectrodeName = [char(i+'A'-1), num2str(j), '-', num2str(m), num2str(n)];
    % recording_duration = max(Times(:)); % Duration spanned by the recording
    clusterIndices = clusterIndices(~cellfun('isempty',clusterIndices));

    % Calculate firing rates and spike counts for each cluster
    clusterStats = struct();

    for clusterIdx = 1:length(clusterIndices)
        templateCount = templateCount + 1;
        validClusterIndices(templateCount) = clusterIdx;
        meanTemplates(:, templateCount) = mean(waveforms(:, clusterIndices{clusterIdx}), 2);
        
        % Calculate statistics for this cluster
        clusterStats(templateCount).numSpikes = length(clusterIndices{clusterIdx});
        
        % Calculate firing rate using actual spike times for this cluster
        t = Times(1, clusterIndices{clusterIdx}); % Get spike start times for this cluster
        t_duration = max(t) - min(t);% Duration spanned by this cluster's spikes
        if length(t) > 1
            
            if t_duration > 0
                clusterStats(templateCount).firingRate = length(clusterIndices{clusterIdx}) / t_duration;
            else
                % If all spikes occur at the same time, use a very small duration
                clusterStats(templateCount).firingRate = length(clusterIndices{clusterIdx}) / 0.001; % 1ms default
            end
        else
            % Single spike case - cannot calculate meaningful firing rate
            clusterStats(templateCount).firingRate = NaN;
        end
        
        clusterStats(templateCount).unitName = [currElectrodeName, '_unit', num2str(clusterIdx)];
        clusterStats(templateCount).duration = t_duration;
        clusterStats(templateCount).firstSpike = min(t);
        clusterStats(templateCount).lastSpike = max(t);
    end
    
    % Normalize templates for better comparison
    templateMin = min(meanTemplates(:));
    templateMax = max(meanTemplates(:));
    normalizedTemplates = (meanTemplates - templateMin) / (templateMax - templateMin);
    
    % Initialize DTW results table
    pairNames = {};
    dtwDistances = [];
    unit1Names = {};
    unit2Names = {};
    unit1NumSpikes = [];
    unit2NumSpikes = [];
    unit1FiringRates = [];
    unit2FiringRates = [];
    unit1Durations = [];
    unit2Durations = [];
    mergeDecisions = {};
    
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
            var_name_unit = [currElectrodeName,'_unit',num2str(template1_idx),'_vs_unit',num2str(template2_idx)];
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
            shapeSimilarityThreshold = 0.3;
            
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
                mergeDecision = 'YES';
            else
                nonMergeCount = nonMergeCount + 1;
                nonMergeDistances(nonMergeCount) = currentDistance;
                mergeDecision = 'NO';
            end
            
            % Create and save figure with enhanced information
            F_electrode = figure('visible','off', 'Position', [100, 100, 800, 600]); 
            
            % Create subplot layout
            subplot(2,1,1);
            hold on;
            plot(waveforms(:,clusterIndices{template1_idx}),'r','LineWidth',0.5, 'Color', [1, 0.7, 0.7]);
            plot(waveforms(:,clusterIndices{template2_idx}),'b','LineWidth',0.5, 'Color', [0.7, 0.7, 1]);
            plot(meanTemplates(:,template1_idx),'r','LineWidth',3);
            plot(meanTemplates(:,template2_idx),'b','LineWidth',3);
            axis([0 40 -1.5e-4 1.5e-4]);
            
            % Enhanced title with DTW distance
            title_text = sprintf('%s\nDTW Distance: %.4f (Threshold: %.3f)', ...
                strrep(var_name_unit, '_', ' '), currentDistance, mergeThreshold);
            title(title_text, 'FontSize', 12, 'FontWeight', 'bold');
            
            legend('Unit 1 spikes', 'Unit 2 spikes', 'Unit 1 mean', 'Unit 2 mean', 'Location', 'best');
            xlabel('Time points');
            ylabel('Amplitude (V)');
            grid on;
            hold off;
            
            % Add statistics subplot
            subplot(2,1,2);
            axis off;
            
            % Prepare statistics text
            unit1_stats = clusterStats(template1_idx);
            unit2_stats = clusterStats(template2_idx);
            
            % Format firing rates (handle NaN values)
            if isnan(unit1_stats.firingRate)
                fr1_text = 'N/A (single spike)';
            else
                fr1_text = sprintf('%.2f Hz', unit1_stats.firingRate);
            end
            
            if isnan(unit2_stats.firingRate)
                fr2_text = 'N/A (single spike)';
            else
                fr2_text = sprintf('%.2f Hz', unit2_stats.firingRate);
            end
            
            % Format peak and trough match results
            if peakMatch
                peakMatchText = 'Yes';
            else
                peakMatchText = 'No';
            end
            
            if troughMatch
                troughMatchText = 'Yes';
            else
                troughMatchText = 'No';
            end
            
            % Create statistics table
            stats_text = {
                '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━';
                sprintf('UNIT COMPARISON STATISTICS');
                '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━';
                '';
                sprintf('%-25s │ %-15s │ %-15s', 'Metric', 'Unit 1 (Red)', 'Unit 2 (Blue)');
                '─────────────────────────┼─────────────────┼─────────────────';
                sprintf('%-25s │ %-15s │ %-15s', 'Unit Name', unit1_stats.unitName, unit2_stats.unitName);
                sprintf('%-25s │ %-15d │ %-15d', 'Number of Spikes', unit1_stats.numSpikes, unit2_stats.numSpikes);
                sprintf('%-25s │ %-15s │ %-15s', 'Firing Rate', fr1_text, fr2_text);
                '';
                sprintf('%-25s │ %-15.3f │ %-15.3f', 'Duration (sec)', unit1_stats.duration, unit2_stats.duration);
                sprintf('%-25s │ %-15.3f │ %-15.3f', 'First Spike (sec)', unit1_stats.firstSpike, unit2_stats.firstSpike);
                sprintf('%-25s │ %-15.3f │ %-15.3f', 'Last Spike (sec)', unit1_stats.lastSpike, unit2_stats.lastSpike);
                '';
                '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━';
                sprintf('DTW ANALYSIS RESULT');
                '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━';
                '';
                sprintf('DTW Distance:     %.4f', currentDistance);
                sprintf('Merge Threshold:  %.4f', mergeThreshold);
                sprintf('Merge Decision:   %s', mergeDecision);
                '';
                sprintf('Peak Match:       %s', peakMatchText);
                sprintf('Trough Match:     %s', troughMatchText);
            };
            
            % Display statistics text
            text(0.05, 0.95, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'FontName', 'FixedWidth', 'FontSize', 9, 'FontWeight', 'normal');
            
            % Add color-coded merge decision box
            if strcmp(mergeDecision, 'YES')
                rectangle('Position', [0.75, 0.1, 0.2, 0.15], 'FaceColor', [0.8, 1, 0.8], ...
                    'EdgeColor', [0, 0.8, 0], 'LineWidth', 2);
                text(0.85, 0.175, {'MERGE', 'YES'}, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0, 0.6, 0]);
            else
                rectangle('Position', [0.75, 0.1, 0.2, 0.15], 'FaceColor', [1, 0.8, 0.8], ...
                    'EdgeColor', [0.8, 0, 0], 'LineWidth', 2);
                text(0.85, 0.175, {'MERGE', 'NO'}, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold', 'FontSize', 12, 'Color', [0.8, 0, 0]);
            end
            
            % Save figure
            saveas(F_electrode,[electrode_DTW_folder,var_name_unit],'png');
            close(F_electrode);
            
            % Store data for Excel table
            pairNames{end+1} = var_name_unit;
            dtwDistances(end+1) = currentDistance;
            unit1Names{end+1} = clusterStats(template1_idx).unitName;
            unit2Names{end+1} = clusterStats(template2_idx).unitName;
            unit1NumSpikes(end+1) = clusterStats(template1_idx).numSpikes;
            unit2NumSpikes(end+1) = clusterStats(template2_idx).numSpikes;
            unit1FiringRates(end+1) = clusterStats(template1_idx).firingRate;
            unit2FiringRates(end+1) = clusterStats(template2_idx).firingRate;
            unit1Durations(end+1) = clusterStats(template1_idx).duration;
            unit2Durations(end+1) = clusterStats(template2_idx).duration;
            mergeDecisions{end+1} = mergeDecision;
        end
    end
    
    % Create DTW results table
    if ~isempty(pairNames)
        T_electrode_DTW = table(pairNames', unit1Names', unit2Names', dtwDistances', ...
            unit1NumSpikes', unit2NumSpikes', unit1FiringRates', unit2FiringRates', ...
            unit1Durations', unit2Durations', mergeDecisions', 'VariableNames', ...
            {'PairName', 'Unit1_Name', 'Unit2_Name', 'DTW_Distance', ...
            'Unit1_NumSpikes', 'Unit2_NumSpikes', 'Unit1_FiringRate_Hz', ...
            'Unit2_FiringRate_Hz', 'Unit1_Duration_sec', 'Unit2_Duration_sec', 'Merge_Decision'});
    else
        % Create empty table with proper variable names
        T_electrode_DTW = table({}, {}, {}, [], [], [], [], [], [], [], {}, 'VariableNames', ...
            {'PairName', 'Unit1_Name', 'Unit2_Name', 'DTW_Distance', ...
            'Unit1_NumSpikes', 'Unit2_NumSpikes', 'Unit1_FiringRate_Hz', ...
            'Unit2_FiringRate_Hz', 'Unit1_Duration_sec', 'Unit2_Duration_sec', 'Merge_Decision'});
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