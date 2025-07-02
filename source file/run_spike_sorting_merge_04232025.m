function run_spike_sorting_merge_04232025(file_folder, threshold_to_merge, refractoryT, ...
    network_participation_threshold, min_spikes_E, max_ISI_E, min_spikes_N, max_ISI_N, ...
    cutoff_frequency, flag_threshold, std_cutoff)
% RUN_SPIKE_SORTING_MERGE_04232025 - Performs spike sorting and burst analysis on neural recordings
%
% This function processes spike data files (.spk), performs clustering, 
% template matching, and analyzes single-unit and network burst activity.
%
% INPUTS:
%   file_folder - Path to folder containing .spk files
%   threshold_to_merge - Distance threshold for merging similar templates (numeric string)
%   refractoryT - Refractory period in milliseconds (numeric string)
%   network_participation_threshold - Threshold for network participation (0-100)
%   min_spikes_E - Minimum number of spikes for electrode burst detection
%   max_ISI_E - Maximum inter-spike interval for electrode burst (ms)
%   min_spikes_N - Minimum number of spikes for network burst detection
%   max_ISI_N - Maximum inter-spike interval for network burst (ms)
%   cutoff_frequency - Firing rate threshold for active vs. inactive units (Hz)
%   flag_threshold - Threshold for flagging potential issues with units
%   std_cutoff - Standard deviation cutoff for outlier removal
%
% OUTPUTS:
%   Creates Excel files and PowerPoint presentations in a 'Results' folder
%   with spike sorting results and burst analysis

    % Initialize parameters
    params = initialize_parameters(threshold_to_merge, refractoryT, ...
        network_participation_threshold, min_spikes_E, max_ISI_E, ...
        min_spikes_N, max_ISI_N, cutoff_frequency, flag_threshold, std_cutoff);
    
    % Find all spike files in the folder
    files = dir(fullfile(file_folder, '*.spk'));
    if isempty(files)
        error('No .spk files found in the specified folder');
    end
    
    % Create output directory
    parentFolder = 'Results';
    if ~exist(parentFolder, 'dir')
        mkdir(parentFolder);
    end
    
    
    
    % Initialize progress bar
    try
        progressbar('Number of Files', 'Progress of current file');
    catch
        warning('Progressbar function not available. Progress will not be displayed.');
    end
    % Save initial variables to avoid clearing them in the loop
    initialVars = who;
    % Process each file
    for fileIdx = 1:length(files)
        if fileIdx == 1
            initialVars{end+1} = 'fileIdx';
        end
        
        % Clear variables but keep initial ones
        clearvars('-except', initialVars{:});
        initialVars = who;
        
        % Process current file
        currentFile = files(fileIdx);
        process_file(currentFile, file_folder, parentFolder, params, fileIdx, length(files));
    end
    
    % Close all figures
    close all;
end




