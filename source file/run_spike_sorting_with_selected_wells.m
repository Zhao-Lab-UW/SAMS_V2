function run_spike_sorting_with_selected_wells(file_folder, selectedWells, ...
    threshold_to_merge, refractoryT, ...
    network_participation_threshold, min_spikes_E, max_ISI_E, ...
    min_spikes_N, max_ISI_N, cutoff_frequency, flag_threshold, std_cutoff)  
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
    
    % Convert selected wells to row/column indices
    wellIndices = cell(length(selectedWells), 4);
    for w = 1:length(selectedWells)
        wellName = selectedWells{w}
        row = double(wellName(1)) - double('A') + 1;
        col = str2double(wellName(2:end));
        wellIndices{w, 1} = row;
        wellIndices{w, 2} = col;
        wellIndices{w, 3} = 1:3; % All electrodes in the well
        wellIndices{w, 4} = 1:3; % All electrodes in the well
    end

    % Process each file
    progressbar('Number of Files', 'Progress of current file');
    for fileIdx = 1:length(files)
        currentFile = files(fileIdx);
        process_file_with_selected_wells_log(currentFile, file_folder, parentFolder, ...
            params, fileIdx, length(files), wellIndices);
    end
    
    close all;
end