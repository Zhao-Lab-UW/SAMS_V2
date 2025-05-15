function [new_idx_list, unit_stats] = process_final_clusters(initial_idx_list, spikes, Times, params, maxTime)
% PROCESS_FINAL_CLUSTERS - Process final clusters after merging and perform burst analysis
%
% INPUTS:
%   initial_idx_list - Cell array of indices for each cluster
%   spikes - Spike waveforms matrix (timepoints × spikes)
%   Times - Spike timing vector
%   params - Parameter structure
%   maxTime - Maximum recording time
%
% OUTPUTS:
%   new_idx_list - Cell array of indices for active units
%   unit_stats - Structure with statistics for each unit

    unit_count = 0;
    new_idx_list = {};
    unit_stats = struct();
    unit_stats.active_units = [];
    unit_stats.inactive_units = [];
    
    % Initialize arrays for firing rates
    unit_stats.active_units_fr = [];
    unit_stats.inactive_units_fr = [];
    
    for k_list = 1:length(initial_idx_list)
        if ~isempty(initial_idx_list{k_list})
            % Get indices for this unit
            index_of_idx = sort(initial_idx_list{k_list});
            
            % Remove outliers based on RMSE
            A = spikes(:, index_of_idx);
            for num = 1:size(A, 1)
                A_RMSE = calculate_WF_RMSE(A);
                [~, TFoutlier] = rmoutliers(A_RMSE, 'mean', 'ThresholdFactor', params.std_cutoff);
                locs = find(TFoutlier == 1);
                A(:, locs) = [];
                index_of_idx(locs) = [];
            end
            
            % Remove spikes with short ISIs (violating refractory period)
            t = Times(1, index_of_idx);
            ISI_s = diff(t);
            burst_pairs_idx = find(ISI_s < params.refractoryT);
            index_of_idx(burst_pairs_idx + 1) = [];
            
            % Update the indices in the initial_idx_list
            initial_idx_list{k_list} = index_of_idx;
            
            % Count number of spikes for this unit
            number_of_spikes = numel(index_of_idx);
            
            % Check if firing rate is above cutoff frequency
            if number_of_spikes / maxTime > params.cutoff_frequency
                % This is an active unit
                unit_count = unit_count + 1;
                
                % Remove invalid time points
                t = Times(1, index_of_idx);
                timeWrong = find(t <= 0);
                index_of_idx(timeWrong) = [];
                
                % Add to active units list
                new_idx_list{unit_count} = index_of_idx;
                
                % Calculate firing rate
                mean_firing_rate = number_of_spikes / (Times(1, max(index_of_idx)) - Times(1, min(index_of_idx)));
                
                % Calculate ISI in general
                ISI_in_general = number_of_spikes / t(end);
                
                % Calculate amplitude statistics
                spikes_amplitude = spikes(:, index_of_idx);
                spikes_amplitude_mean = mean(spikes_amplitude, 2);
                [Amplitude_Avg, Amplitude_Avg_loc] = max(abs(spikes_amplitude_mean));
                Amplitude_Avg = Amplitude_Avg * 1000; % Convert to mV
                Amplitude_Std = std(spikes_amplitude(Amplitude_Avg_loc, :)) * 1000;
                
                % Calculate peak-to-trough amplitude
                [spikes_amplitude_mean_peak, spikes_amplitude_mean_peak_locs] = max(spikes_amplitude_mean);
                [spikes_amplitude_mean_trough, spikes_amplitude_mean_trough_locs] = min(spikes_amplitude_mean);
                Amplitude_Avg_Peak_to_Trough = (spikes_amplitude_mean_peak - spikes_amplitude_mean_trough) * 1000;
                Amplitude_Std_Peak_to_Trough = std(spikes_amplitude(spikes_amplitude_mean_peak_locs, :) - ...
                    spikes_amplitude(spikes_amplitude_mean_trough_locs, :)) * 1000;
                
                % Calculate ISI in ms
                ISI_ms = diff(t) * 1000;
                shortest_ISI_ms = min(ISI_ms);
                
                % Get burst information
                burst_info = get_burst(t, params.maximum_inter_spike_interval, params.minimum_spike_per_burst);
                
                % Store unit statistics
                if isempty(burst_info)
                    % No bursts detected
                    unit_stats.active_units(unit_count).stats = [number_of_spikes; mean_firing_rate; ...
                        ISI_in_general; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; ...
                        Amplitude_Avg; Amplitude_Std; Amplitude_Avg_Peak_to_Trough; ...
                        Amplitude_Std_Peak_to_Trough; shortest_ISI_ms; length(burst_pairs_idx)];
                else
                    % Bursts detected, calculate burst statistics
                    burst_duration_avg = mean(burst_info(3, :));
                    burst_duration_std = std(burst_info(3, :));
                    burst_duration_median = median(burst_info(3, :));
                    normalized_duration_IQR = iqr(burst_info(3, :)) / burst_duration_median;
                    burst_counts = size(burst_info, 2);
                    burst_percentage = sum(burst_info(2, :)) / number_of_spikes * 100;
                    burst_frequency = burst_counts / t(end);
                    number_of_spikes_per_burst_avg = mean(burst_info(2, :));
                    number_of_spikes_per_burst_std = std(burst_info(2, :));
                    number_of_spikes_per_burst_median = median(burst_info(2, :));
                    mean_ISI_within_burst_avg = mean(burst_info(4, :));
                    mean_ISI_within_burst_std = std(burst_info(4, :));
                    mean_ISI_within_burst_median = median(burst_info(4, :));
                    ISI_CV = mean_ISI_within_burst_std / mean_ISI_within_burst_avg;
                    
                    % Calculate inter-burst interval statistics
                    inter_burst_interval_avg = mean(diff(t(burst_info(1, :))));
                    inter_burst_interval_std = std(diff(t(burst_info(1, :))));
                    inter_burst_interval_median = median(diff(t(burst_info(1, :))));
                    inter_burst_interval_CV = inter_burst_interval_std / inter_burst_interval_avg;
                    
                    % Store stats in the same order as the table expects
                    unit_stats.active_units(unit_count).stats = [number_of_spikes; mean_firing_rate; ...
                        ISI_in_general; burst_counts; burst_frequency; burst_duration_avg; ...
                        burst_duration_std; burst_duration_median; burst_percentage; ...
                        normalized_duration_IQR; number_of_spikes_per_burst_avg; ...
                        number_of_spikes_per_burst_std; number_of_spikes_per_burst_median; ...
                        mean_ISI_within_burst_avg; mean_ISI_within_burst_std; ...
                        mean_ISI_within_burst_median; ISI_CV; inter_burst_interval_avg; ...
                        inter_burst_interval_std; inter_burst_interval_median; ...
                        inter_burst_interval_CV; Amplitude_Avg; Amplitude_Std; ...
                        Amplitude_Avg_Peak_to_Trough; Amplitude_Std_Peak_to_Trough; ...
                        shortest_ISI_ms; length(burst_pairs_idx)];
                end
                
                % Store the unit index in the stats
                unit_stats.active_units(unit_count).index = k_list;
                
                % Add to active units firing rate list
                unit_stats.active_units_fr(unit_count) = number_of_spikes / maxTime;
            else
                % This is an inactive unit (firing rate below cutoff)
                inactive_idx = length(unit_stats.inactive_units) + 1;
                unit_stats.inactive_units(inactive_idx).index = k_list;
                unit_stats.inactive_units_fr(inactive_idx) = number_of_spikes / maxTime;
            end
        end
    end
end

% function [new_idx_list, unit_stats] = process_final_clusters(initial_idx_list, spikes, Times, params, maxTime)
%     % Process final clusters after merging and perform burst analysis
%     %
%     % INPUTS:
%     %   initial_idx_list - Cell array of indices for each cluster
%     %   spikes - Spike waveforms matrix (time points × spikes)
%     %   Times - Spike timing vector
%     %   params - Parameter structure
%     %   maxTime - Maximum recording time
%     %
%     % OUTPUTS:
%     %   new_idx_list - Cell array of indices for active units
%     %   unit_stats - Structure with statistics for each unit
% 
%     unit_count = 0;
%     new_idx_list = {};
%     unit_stats = struct();
%     unit_stats.active_units = [];
%     unit_stats.inactive_units = [];
% 
%     for k_list = 1:length(initial_idx_list)
%         if ~isempty(initial_idx_list{k_list})
%             % Get indices for this unit
%             index_of_idx = sort(initial_idx_list{k_list});
% 
%             % Remove outliers based on RMSE
%             A = spikes(:, index_of_idx);
%             for num = 1:size(A, 1)
%                 A_RMSE = calculate_WF_RMSE(A);
%                 [~, TFoutlier] = rmoutliers(A_RMSE, 'mean', 'ThresholdFactor', params.std_cutoff);
%                 locs = find(TFoutlier == 1);
%                 A(:, locs) = [];
%                 index_of_idx(locs) = [];
%             end
% 
%             % Remove spikes with short ISIs (violating refractory period)
%             t = Times(1, index_of_idx);
%             ISI_s = diff(t);
%             burst_pairs_idx = find(ISI_s < params.refractoryT);
%             index_of_idx(burst_pairs_idx + 1) = [];
% 
%             % Update the indices in the initial_idx_list
%             initial_idx_list{k_list} = index_of_idx;
% 
%             % Count number of spikes for this unit
%             number_of_spikes = numel(index_of_idx);
% 
%             % Check if firing rate is above cutoff frequency
%             if number_of_spikes / maxTime > params.cutoff_frequency
%                 % This is an active unit
%                 unit_count = unit_count + 1;
% 
%                 % Remove invalid time points
%                 t = Times(1, index_of_idx);
%                 timeWrong = find(t <= 0);
%                 index_of_idx(timeWrong) = [];
% 
%                 % Add to active units list
%                 new_idx_list{unit_count} = index_of_idx;
% 
%                 % Calculate firing rate
%                 mean_firing_rate = number_of_spikes / (Times(1, max(index_of_idx)) - Times(1, min(index_of_idx)));
% 
%                 % Calculate ISI in general
%                 ISI_in_general = number_of_spikes / t(end);
% 
%                 % Calculate amplitude statistics
%                 spikes_amplitude = spikes(:, index_of_idx);
%                 spikes_amplitude_mean = mean(spikes_amplitude, 2);
%                 [Amplitude_Avg, Amplitude_Avg_loc] = max(abs(spikes_amplitude_mean));
%                 Amplitude_Avg = Amplitude_Avg * 1000; % Convert to mV
%                 Amplitude_Std = std(spikes_amplitude(Amplitude_Avg_loc, :)) * 1000;
% 
%                 % Calculate peak-to-trough amplitude
%                 [spikes_amplitude_mean_peak, spikes_amplitude_mean_peak_locs] = max(spikes_amplitude_mean);
%                 [spikes_amplitude_mean_trough, spikes_amplitude_mean_trough_locs] = min(spikes_amplitude_mean);
%                 Amplitude_Avg_Peak_to_Trough = (spikes_amplitude_mean_peak - spikes_amplitude_mean_trough) * 1000;
%                 Amplitude_Std_Peak_to_Trough = std(spikes_amplitude(spikes_amplitude_mean_peak_locs, :) - ...
%                     spikes_amplitude(spikes_amplitude_mean_trough_locs, :)) * 1000;
% 
%                 % Calculate ISI in ms
%                 ISI_ms = diff(t) * 1000;
%                 shortest_ISI_ms = min(ISI_ms);
% 
%                 % Get burst information
%                 burst_info = get_burst(t, params.maximum_inter_spike_interval, params.minimum_spike_per_burst);
% 
%                 % Store unit statistics
%                 if isempty(burst_info)
%                     % No bursts detected
%                     unit_stats.active_units(unit_count).stats = [number_of_spikes; mean_firing_rate; ...
%                         ISI_in_general; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; ...
%                         Amplitude_Avg; Amplitude_Std; Amplitude_Avg_Peak_to_Trough; ...
%                         Amplitude_Std_Peak_to_Trough; shortest_ISI_ms; length(burst_pairs_idx)];
%                 else
%                     % Bursts detected, calculate burst statistics
%                     burst_duration_avg = mean(burst_info(3, :));
%                     burst_duration_std = std(burst_info(3, :));
%                     burst_duration_median = median(burst_info(3, :));
%                     normalized_duration_IQR = iqr(burst_info(3, :)) / burst_duration_median;
%                     burst_counts = size(burst_info, 2);
%                     burst_percentage = sum(burst_info(2, :)) / number_of_spikes * 100;
%                     burst_frequency = burst_counts / t(end);
%                     number_of_spikes_per_burst_avg = mean(burst_info(2, :));
%                     number_of_spikes_per_burst_std = std(burst_info(2, :));
%                     number_of_spikes_per_burst_median = median(burst_info(2, :));
%                     mean_ISI_within_burst_avg = mean(burst_info(4, :));
%                     mean_ISI_within_burst_std = std(burst_info(4, :));
%                     mean_ISI_within_burst_median = median(burst_info(4, :));
%                     ISI_CV = mean_ISI_within_burst_std / mean_ISI_within_burst_avg;
% 
%                     % Calculate inter-burst interval statistics
%                     inter_burst_interval_avg = mean(diff(t(burst_info(1, :))));
%                     inter_burst_interval_std = std(diff(t(burst_info(1, :))));
%                     inter_burst_interval_median = median(diff(t(burst_info(1, :))));
%                     inter_burst_interval_CV = inter_burst_interval_std / inter_burst_interval_avg;
% 
%                     % Store stats in the same order as the table expects
%                     unit_stats.active_units(unit_count).stats = [number_of_spikes; mean_firing_rate; ...
%                         ISI_in_general; burst_counts; burst_frequency; burst_duration_avg; ...
%                         burst_duration_std; burst_duration_median; burst_percentage; ...
%                         normalized_duration_IQR; number_of_spikes_per_burst_avg; ...
%                         number_of_spikes_per_burst_std; number_of_spikes_per_burst_median; ...
%                         mean_ISI_within_burst_avg; mean_ISI_within_burst_std; ...
%                         mean_ISI_within_burst_median; ISI_CV; inter_burst_interval_avg; ...
%                         inter_burst_interval_std; inter_burst_interval_median; ...
%                         inter_burst_interval_CV; Amplitude_Avg; Amplitude_Std; ...
%                         Amplitude_Avg_Peak_to_Trough; Amplitude_Std_Peak_to_Trough; ...
%                         shortest_ISI_ms; length(burst_pairs_idx)];
%                 end
% 
%                 % Store the unit index in the stats
%                 unit_stats.active_units(unit_count).index = k_list;
% 
%                 % Add to active units firing rate list
%                 unit_stats.active_units_fr(unit_count) = number_of_spikes / maxTime;
%             else
%                 % This is an inactive unit (firing rate below cutoff)
%                 unit_stats.inactive_units(end+1).index = k_list;
%                 unit_stats.inactive_units_fr(end+1) = number_of_spikes / maxTime;
%             end
%         end
%     end
% end