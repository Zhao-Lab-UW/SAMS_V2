function flagTooMuch = check_too_many_outliers(initial_idx_list, idx, flag_threshold)
    % Check if too many outliers were removed from a cluster
    %
    % INPUTS:
    %   initial_idx_list - Cell array of indices for each cluster after outlier removal
    %   idx - Original cluster indices from spectral clustering
    %   flag_threshold - Threshold ratio for flagging excessive outlier removal
    %
    % OUTPUTS:
    %   flagTooMuch - Boolean flag indicating if too many outliers were removed

    flagTooMuch = false;
    
    for k_file = 1:length(initial_idx_list)
        if ~isempty(initial_idx_list{k_file})
            % Count original spikes in this cluster
            orig_count = sum(idx == k_file);
            
            % Count spikes after outlier removal
            new_count = length(initial_idx_list{k_file});
            
            % Calculate the fraction of removed spikes
            removed_fraction = 1 - (new_count / orig_count);
            
            % Flag if too many spikes were removed
            if removed_fraction >= flag_threshold
                flagTooMuch = true;
                return;
            end
        end
    end
end