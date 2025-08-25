function initial_idx_list = process_initial_clusters(idx, spikes, std_cutoff)
    % Process initial clusters from spectral clustering and remove outliers
    %
    % INPUTS:
    %   idx - Cluster indices from spectral clustering
    %   spikes - Spike waveforms matrix (time points × spikes)
    %   std_cutoff - Standard deviation threshold for outlier removal
    %
    % OUTPUTS:
    %   initial_idx_list - Cell array of indices for each cluster after outlier removal

    initial_idx_list = cell(1, max(idx));
    outlier_WFs = [];

    for k_file = 1:max(idx)
        % Get indices of spikes in this cluster
        index_of_idx = find(idx == k_file);
        
        % Get spike waveforms for this cluster
        A = spikes(:, index_of_idx);
        
        % Track number of spikes before removing outliers
        num_spikes_before_removing_outliers = numel(index_of_idx);
        
        % Calculate RMSE for each waveform compared to mean
        A_RMSE = calculate_WF_RMSE(A);
        
        % Remove outliers based on RMSE
        [~, TFoutlier] = rmoutliers(A_RMSE, 'mean', 'ThresholdFactor', std_cutoff);
        locs = find(TFoutlier == 1);
        
        % Remove outlier waveforms
        index_of_idx(locs) = [];
        
        % Store indices of non-outlier spikes for this cluster
        initial_idx_list{k_file} = index_of_idx;
    end
end