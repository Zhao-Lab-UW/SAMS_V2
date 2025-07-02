function [initial_idx_list, HDT_flag] = fix_undersorting(initial_idx_list, spikes, Times, score, refractoryT, params)
    % Fix undersorting using Hartigan's dip test for bimodality
    %
    % INPUTS:
    %   initial_idx_list - Cell array of indices for each cluster after outlier removal
    %   spikes - Spike waveforms matrix (time points × spikes)
    %   Times - Spike timing vector
    %   score - PCA scores for each spike
    %   refractoryT - Refractory period in seconds
    %   params - Parameter structure
    %
    % OUTPUTS:
    %   initial_idx_list - Updated cell array of indices after fixing undersorting
    %   HDT_flag - Flag indicating if Hartigan's dip test found bimodality

    HDT_flag = 0;
    nboot = 500; % Bootstrap sample size for the dip test
    count_fix_undersorting = 0;
    copy_initial_idx_list = initial_idx_list;
    initial_idx_list = cell(1, length(copy_initial_idx_list)*2); % Pre-allocate larger array
    
    for num_unit_initial = 1:length(copy_initial_idx_list)
        % try
            if length(copy_initial_idx_list{num_unit_initial}) > 2
                % Get PCA scores for this cluster
                samplePCA1 = score(copy_initial_idx_list{num_unit_initial}, 1)';
                samplePCA2 = score(copy_initial_idx_list{num_unit_initial}, 2)';
                
                % Perform Hartigan's dip test on the first two PCs
                [~, p_value1, ~, ~] = HartigansDipSignifTest(samplePCA1, nboot);
                [~, p_value2, ~, ~] = HartigansDipSignifTest(samplePCA2, nboot);

                % If either PC shows significant bimodality (p < 0.05)
                if (p_value1 < 0.05) || (p_value2 < 0.05)
                    try
                        % Perform additional clustering on this unit
                        sample2d = score(copy_initial_idx_list{num_unit_initial}, 1:2);
                        [bimodal_idx, ~] = spectralcluster(sample2d, 2);
                    catch ME
                        if contains(ME.message, 'Invalid data type') || contains(ME.message, 'real array')
                            try
                                % Try k-means as fallback
                                [bimodal_idx, ~] = kmeans(sample2d', 2, 'Replicates', 2);
                                warning('Used k-means fallback');
                            catch
                                % If even k-means fails, skip this iteration
                                warning('Both spectral clustering and k-means failed');
                                continue;
                            end
                        else
                            rethrow(ME);  % Re-throw if it's a different error
                        end
                    end
                    % Split the cluster into two based on the new clustering
                    bimodal_idx1 = copy_initial_idx_list{num_unit_initial}(bimodal_idx == 1);
                    bimodal_idx2 = copy_initial_idx_list{num_unit_initial}(bimodal_idx == 2);

                    % Get spike waveforms and times for each sub-cluster
                    spk1 = spikes(:, bimodal_idx1); 
                    t1 = Times(:, bimodal_idx1);
                    spk2 = spikes(:, bimodal_idx2); 
                    t2 = Times(:, bimodal_idx2);
                    
                    % Calculate overlap between the two sub-clusters
                    [ovlp1, ovlp2, tdiff] = calculate_overlap(spk1, spk2, t1, t2);
                    ovlp = max(sum(ovlp2 > 0.7) / length(ovlp1), sum(ovlp2 > 0.7) / length(ovlp2));
                    
                    % Check if the sub-clusters are likely from the same spike
                    sameSpike = (ovlp > 0.7) && (tdiff(1) > refractoryT);
                    
                    if sameSpike
                        % Keep as a single unit if they're the same spike
                        count_fix_undersorting = count_fix_undersorting + 1;
                        initial_idx_list{count_fix_undersorting} = copy_initial_idx_list{num_unit_initial};
                    else
                        % Split into two units if they're different spikes
                        HDT_flag = 1;
                        count_fix_undersorting = count_fix_undersorting + 1;
                        initial_idx_list{count_fix_undersorting} = bimodal_idx1;
                        count_fix_undersorting = count_fix_undersorting + 1;
                        initial_idx_list{count_fix_undersorting} = bimodal_idx2;
                    end
                else
                    % Keep as a single unit if no bimodality detected
                    count_fix_undersorting = count_fix_undersorting + 1;
                    initial_idx_list{count_fix_undersorting} = copy_initial_idx_list{num_unit_initial};
                end
            end
        % catch
        %     % Keep as is if an error occurs
        %     fprintf('Error in dip test for unit %d, keeping as is\n', num_unit_initial);
        %     count_fix_undersorting = count_fix_undersorting + 1;
        %     initial_idx_list{count_fix_undersorting} = copy_initial_idx_list{num_unit_initial};
        % end
    end
    
    % Trim empty cells
    initial_idx_list = initial_idx_list(1:count_fix_undersorting);
end