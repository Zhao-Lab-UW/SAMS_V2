function [electrode_results, electrode_stats] = process_electrode(spikes, Times, i, j, m, n, params, maxTime, maxSpikeAmp, minSpikeAmp, score)
% PROCESS_ELECTRODE - Process a single electrode
%
% INPUTS:
%   spikes - Spike waveforms matrix (timepoints × spikes)
%   Times - Spike timing vector
%   i, j, m, n - Electrode indices
%   params - Parameter structure
%   maxTime - Maximum recording time
%   maxSpikeAmp, minSpikeAmp - Amplitude range
%   score - PCA scores (optional)
%
% OUTPUTS:
%   electrode_results - Structure with electrode results
%   electrode_stats - Structure with electrode statistics

    % Initialize results
    electrode_results = struct();
    electrode_stats = struct();
    electrode_results.new_idx_list = {};
    electrode_results.HDT_flag = 0;
    electrode_results.DTW_flag = 0;
    electrode_results.possibleMUA = 0;
    electrode_results.flagTooMuch = 0;
    
    % Current electrode name
    currElectrodeName = [char(i+'A'-1), num2str(j), '_', num2str(m), num2str(n)];
    electrode_results.name = currElectrodeName;
    
    % Check if enough spikes
    numSpikes = size(spikes, 2);
    if numSpikes < 5
        fprintf('Too few spikes (%d) for processing.\n', numSpikes);
        return;
    end
    
    % Perform PCA if not provided
    if nargin < 11 || isempty(score)
        try
            [~, score, ~] = pca(spikes');
        catch ME
            fprintf('Error in PCA: %s\n', ME.message);
            return;
        end
    end
    
    % Check if enough dimensions in score
    if size(score, 2) < 2
        fprintf('Not enough dimensions in PCA score.\n');
        return;
    end
    
    % Initialize spectral clustering parameters
    myfunc = @(X, K)(spectralcluster(X, K));
    klist = 1:5; % Number of clusters to try
    
    % Evaluate optimal number of clusters
    try
        eva = evalclusters(score(:,1:2), myfunc, "DaviesBouldin", 'klist', klist);
        
        if isnan(eva.OptimalK)
            fprintf('Could not determine optimal number of clusters. Using K=2.\n');
            K = 2;
        else
            fprintf('Optimal number of clusters: %d\n', eva.OptimalK);
            K = eva.OptimalK;
        end
    catch ME
        fprintf('Error in cluster evaluation: %s\n', ME.message);
        fprintf('Using default K=2.\n');
        K = 2;
    end
    
    % Perform spectral clustering
    try
        if K > 1
            [idx, ~] = spectralcluster(score(:,1:2), K);
        else
            idx = ones(size(score, 1), 1);
        end
    catch ME
        fprintf('Error in spectral clustering: %s\n', ME.message);
        idx = ones(size(score, 1), 1);
    end
    
    % Process initial clusters
    initial_idx_list = process_initial_clusters(idx, spikes, params.std_cutoff);
    
    % Check for too many outliers
    flagTooMuch = check_too_many_outliers(initial_idx_list, idx, params.flag_threshold);
    electrode_results.flagTooMuch = flagTooMuch;
    if flagTooMuch
        electrode_results.flagTooMuchName = currElectrodeName;
    end
    
    % Fix under-sorting
    [initial_idx_list, HDT_flag] = fix_undersorting(initial_idx_list, spikes, Times, score, params.refractoryT, params);
    electrode_results.HDT_flag = HDT_flag;
    
    % Perform template comparison and merge similar units
    if length(initial_idx_list) > 1
        [dist, merge, distMerge, distNotMerge] = template_comparison_01042024(spikes, initial_idx_list, params.threshold_to_merge);
        
        % Process outliers if needed
        if electrode_results.flagTooMuch == 1
            try
                outlier_WFs = setdiff(1:size(spikes, 2), [initial_idx_list{:}]);
                [merge_outlier_list, initial_idx_list_unsorted_WFs] = template_comparison_outlier_02212024(spikes, initial_idx_list, outlier_WFs, params.threshold_to_merge, myfunc, klist, score, params.std_cutoff);
                
                % Add outliers back to clusters
                if ~isempty(merge_outlier_list)
                    for ii = 1:length(initial_idx_list)
                        initial_idx_list{ii} = sort([initial_idx_list{ii}; merge_outlier_list{ii}]);
                    end
                end
            catch ME
                fprintf('Error processing outliers: %s\n', ME.message);
                initial_idx_list_unsorted_WFs = {};
            end
        else
            initial_idx_list_unsorted_WFs = {};
        end
        
        % Merge units based on DTW distance
        if ~isempty(merge)
            electrode_results.DTW_flag = 1;
            for ii = 1:size(merge, 2)
                for jj = 2:size(merge{ii}, 2)
                    initial_idx_list{merge{ii}(1)} = sort([initial_idx_list{merge{ii}(1)}; initial_idx_list{merge{ii}(jj)}]);
                    initial_idx_list{merge{ii}(jj)} = [];
                end
            end
        end
        
        % Add unsorted waveforms
        initial_idx_list = [initial_idx_list, initial_idx_list_unsorted_WFs];
        
        % Perform second merge if needed
        if length(initial_idx_list) > 1
            [dist, merge, distMerge, distNotMerge] = template_comparison_01042024(spikes, initial_idx_list, params.threshold_to_merge);
            
            if ~isempty(merge)
                electrode_results.DTW_flag = 1;
                for ii = 1:size(merge, 2)
                    for jj = 2:size(merge{ii}, 2)
                        initial_idx_list{merge{ii}(1)} = sort([initial_idx_list{merge{ii}(1)}; initial_idx_list{merge{ii}(jj)}]);
                        initial_idx_list{merge{ii}(jj)} = [];
                    end
                end
            end
        end
    end
    
    % Process final clusters
    [new_idx_list, unit_stats] = process_final_clusters(initial_idx_list, spikes, Times, params, maxTime);
    electrode_results.new_idx_list = new_idx_list;
    electrode_stats = unit_stats;
    
    % Check for possible multi-unit activity
    if length(new_idx_list) > 1
        [dist, merge, distMerge, distNotMerge] = template_comparison_01042024(spikes, initial_idx_list, params.threshold_to_merge);
        if ~isempty(find(distNotMerge < params.threshold_to_merge*1.5, 1))
            electrode_results.possibleMUA = 1;
            electrode_results.possibleMUAName = currElectrodeName;
        end
    end
    
    % Check for inconsistency in flags
    if electrode_results.HDT_flag ~= electrode_results.DTW_flag
        fprintf('Flag inconsistency in electrode %s: HDT=%d, DTW=%d\n', currElectrodeName, electrode_results.HDT_flag, electrode_results.DTW_flag);
    end
end