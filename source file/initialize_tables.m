function [T, T_electrode, T_parameters] = initialize_tables(params)
    % Create table for unit measurements
    Measurements = {'Number of Spikes', 'Mean Firing Rate (Hz)', 'ISI in general',...
        'Number of Bursts',...
        'Burst Frequency (Hz)',...
        'Burst Duration -Avg(sec)','Burst Duration -Std(sec)', 'Burst Duration -Median(sec)', ...
        'Burst Percentage (%)',...
        'Normalized Burst Duration IQR',...
        'Number of Spikes per Burst -Avg', 'Number of Spikes per Burst -Std', 'Number of Spikes per Burst -Median',...
        'Mean ISI within Burst -Avg', 'Mean ISI within Burst -Std', 'Mean ISI within Burst -Median',...
        'ISI Coefficient of Variation',...
        'Inter-Burst Interval - Avg', 'Inter-Burst Interval - Std', 'Inter-Burst Interval - Median',...
        'Inter-Burst Interval Coefficient of Variation',...
        'Amplitude -Avg(mV)','Amplitude -Std(mV)',...
        'Amplitude - Avg Peak to Trough(mV)', 'Amplitude - Std Peak to Trough(mV)',...
        'max ISI (ms)','number of removed short ISIs'};
    T = table(Measurements');
    
    % Create expanded table for electrode measurements - fixed array concatenation
    Measurements_electrode = {
        'Total number of spikes'; 
        'Number of excluded spikes';
        'Unsorted firing rate (Hz)'; 
        'Number of Units Detected';
        'Average spike frequency per electrode (Hz)';
        'Electrode firing rate (all units combined, Hz)';
        'Unit firing rates - Mean (Hz)';
        'Unit firing rates - Std (Hz)';
        'Unit firing rates - Median (Hz)';
        'Unit firing rates - Min (Hz)';
        'Unit firing rates - Max (Hz)';
        'Unit spike counts - Mean';
        'Unit spike counts - Std';
        'Unit spike counts - Median';
        'Unit spike counts - Min';
        'Unit spike counts - Max';
        'Number of bursting units';
        'Percentage of bursting units (%)';
        'Average bursts per unit';
        'Average burst duration per electrode (sec)';
        'Average spikes per burst per electrode';
        'Inter-burst interval - electrode average (sec)';
        'Electrode burst frequency (Hz)';
        'Total recording duration for electrode (sec)';
        'Electrode activity span (first to last spike, sec)';
        'Electrode activity percentage (%)'
    };
    T_electrode = table(Measurements_electrode);
    
    % Create parameters table
    parameter_list = {'threshold_to_merge', 'refractoryT', ...
        'min_spikes_E', 'max_ISI_E', 'min_spikes_N', 'max_ISI_N',...
        'network_participation_threshold',...
        'cutoff_frequency', 'flag_threshold', 'std_cutoff'};
    T_parameters = table(parameter_list');
    
    % Add parameter values
    parameter_values = [params.threshold_to_merge; params.refractoryT; ...
        params.min_spikes_E; params.max_ISI_E; params.min_spikes_N; params.max_ISI_N; ...
        params.network_participation_threshold; ...
        params.cutoff_frequency; params.flag_threshold; params.std_cutoff];
    Tleft_parameters = table(parameter_values);
    Tleft_parameters.Properties.VariableNames(1) = {'parameter Value'};
    T_parameters = [T_parameters, Tleft_parameters];
end