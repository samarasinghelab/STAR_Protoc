%% Calcium Imaging Analysis Script
% This script analyzes processed calcium imaging data from CaImAn/ezcalcium
% It calculates burstiness, firing rates, spike probabilities, functional connectivity, and synchronization metrics

%% Clear workspace
clear;
clc;

%% User-defined parameters (modify these as needed)
Nsamples = 1000;         % Number of samples for spike train generation
frame_rate = 0.3875;     % Frame rate in seconds per frame
max_lag = 13;            % Maximum lag for cross-correlation in frames (~5 seconds at default frame rate)

% Statistical significance parameters
n_shuffles = 1000;       % Number of shuffles for significance testing
significance_alpha = 0.10; % Significance level (e.g., 0.10 for 90% confidence)

%% Select files to process
fprintf('\nSelect processed files to analyze...\n');
[filenames, folderpath] = uigetfile('*_mcor.mat', ...
                                     'Select _mcor.mat files to analyze', ...
                                     'MultiSelect', 'on');
if isequal(filenames, 0)
    error('No files selected. Exiting.');
end

% Handle single vs multiple file selection
if ischar(filenames)
    filenames = {filenames};  % Convert to cell array for consistency
end

filelist = fullfile(folderpath, filenames);
fprintf('Selected %d files for processing.\n\n', numel(filelist));

%% Process each file
for fi = 1:numel(filelist)
    fil = filelist{fi};
    fprintf('Processing file %d/%d: %s\n', fi, numel(filelist), fil);
    
    %% Load data
    % Load main variables
    lf = load(fil, 'ROI_center_refined', 'F_raw_refined', ...
              'F_inferred_refined', 'Z_mod_refined');
    
    % Try to load parameters used for Z_mod_refined calculation
    params = struct();
    try
        lf_params = load(fil, 'params');
        params = lf_params.params;
        fprintf('  Loaded existing parameters from file.\n');
    catch
        fprintf('  No parameters found in file. Using defaults.\n');
    end
    
    % Add Nsamples to parameters
    params.Nsamples = Nsamples;
    
    [n_rois, n_frames] = size(lf.Z_mod_refined);
    total_time = n_frames * frame_rate; % Total recording time in seconds
    
    fprintf('  Number of ROIs: %d\n', n_rois);
    fprintf('  Number of frames: %d\n', n_frames);
    fprintf('  Total time: %.2f seconds\n', total_time);
    
    %% Initialize output variables
    burstiness = zeros(n_rois, 1);
    firing_rate_hz = zeros(n_rois, 1);
    P_spike = zeros(n_rois, n_frames);
    
    %% Process each ROI
    fprintf('  Analyzing ROI spike statistics...\n');
    for i = 1:n_rois
        % Generate N spike trains for this ROI
        sinfo = cont_ca_sampler(lf.Z_mod_refined(i, :), params);
        ntrials = size(sinfo.ss, 1);
        
        % Collect all spike times across trials
        all_isi = [];  % Inter-spike intervals
        total_spikes = 0;
        
        % Calculate spike probability (histogram across bins)
        spike_counts = zeros(1, n_frames);
        
        for trial = 1:ntrials
            spike_times = sinfo.ss{trial};
            
            % Count spikes in each frame bin
            spike_counts = spike_counts + histcounts(spike_times, 0:n_frames);
            
            % Calculate inter-spike intervals for this trial
            if numel(spike_times) > 1
                isi = diff(spike_times) * frame_rate; % Convert to seconds
                all_isi = [all_isi; isi(:)];
            end
            
            total_spikes = total_spikes + numel(spike_times);
        end
        
        % Calculate burstiness (coefficient of variation of ISI)
        if ~isempty(all_isi) && std(all_isi) > 0
            burstiness(i) = std(all_isi) / mean(all_isi);
        else
            burstiness(i) = 0;
        end
        
        % Calculate firing rate in Hz
        firing_rate_hz(i) = total_spikes / (ntrials * total_time);
        
        % Calculate spike probability (average across trials)
        P_spike(i, :) = spike_counts / ntrials;
    end
    
    fprintf('  Burstiness range: [%.3f, %.3f]\n', min(burstiness), max(burstiness));
    fprintf('  Firing rate range: [%.3f, %.3f] Hz\n', min(firing_rate_hz), max(firing_rate_hz));
    
    %% Calculate functional connectivity using cross-correlation
    fprintf('  Computing functional connectivity matrix...\n');
    
    % Initialize matrices
    R = zeros(n_rois, n_rois);              % Cross-correlation values
    R_pvalue = ones(n_rois, n_rois);        % P-values for significance
    R_significant = false(n_rois, n_rois); % Boolean mask for significant connections
    
    % Compute null distribution using spike train shuffling
    fprintf('  Computing significance thresholds via spike train shuffling (%d shuffles)...\n', n_shuffles);
    
    for i = 1:n_rois
        for j = i:n_rois  % Only compute upper triangle (symmetric)
            if i == j
                R(i, j) = 1.0;  % Perfect autocorrelation
                R_pvalue(i, j) = 0;
                R_significant(i, j) = true;
            else
                % Compute actual cross-correlation
                [xcorr_vals, lags] = xcorr(P_spike(i, :), P_spike(j, :), max_lag, 'coeff');
                
                % Store signed maximum correlation value (keep sign for excitatory/inhibitory)
                [~, max_idx] = max(abs(xcorr_vals));
                actual_corr = xcorr_vals(max_idx);
                R(i, j) = actual_corr;
                R(j, i) = actual_corr;  % Symmetric
                
                % Generate null distribution by shuffling one spike train
                null_corr = zeros(n_shuffles, 1);
                for shuffle = 1:n_shuffles
                    % Randomly circularly shift one spike train to destroy temporal relationship
                    random_shift = randi([1, n_frames - 1]);
                    P_shuffled = circshift(P_spike(j, :), random_shift);
                    
                    % Compute cross-correlation on shuffled data
                    xcorr_jitter = xcorr(P_spike(i, :), P_shuffled, max_lag, 'coeff');
                    [~, max_idx_jitter] = max(abs(xcorr_jitter));
                    null_corr(shuffle) = xcorr_jitter(max_idx_jitter);
                end
                
                % Calculate p-value (two-tailed test)
                p_val = sum(abs(null_corr) >= abs(actual_corr)) / n_shuffles;
                R_pvalue(i, j) = p_val;
                R_pvalue(j, i) = p_val;  % Symmetric
                
                % Determine significance
                is_significant = p_val < significance_alpha;
                R_significant(i, j) = is_significant;
                R_significant(j, i) = is_significant;  % Symmetric
            end
        end
    end
    
    %% Calculate synchronization metrics
    % Count significant connections
    upper_triangle_indices = triu(true(n_rois), 1);  % Exclude diagonal
    unique_correlations = R(upper_triangle_indices);
    n_significant = sum(R_significant(upper_triangle_indices));
    percent_significant = 100 * n_significant / numel(unique_correlations);
    
    % Global population synchrony (multivariate measure)
    % Uses ALL neurons - measures population-level coordination regardless of pairwise significance
    fprintf('  Computing global population synchrony...\n');
    % Normalize each neuron (z-score) so all neurons contribute equally
    P_spike_norm = zscore(P_spike, 0, 2);  % Normalize each row (neuron) across time
    % At each time point, measure standard deviation across normalized neurons
    global_synchrony_timeseries = 1 ./ (std(P_spike_norm, 0, 1) + eps);  % Inverse std at each time
    % Summary index: average synchrony across time
    global_synchrony_index = mean(global_synchrony_timeseries);
    
    % Spectral coherence (frequency-domain synchronization)
    % Computes for ALL pairs - can filter by R_significant later for analysis
    fprintf('  Computing spectral coherence...\n');
    coherence_matrix = zeros(n_rois, n_rois);
    for i = 1:n_rois
        coherence_matrix(i, i) = 1.0;  % Perfect coherence with self
        for j = i+1:n_rois
            % Magnitude-squared coherence
            [Cxy, ~] = mscohere(P_spike(i, :), P_spike(j, :), [], [], [], 1/frame_rate);
            % Average coherence across frequencies
            avg_coherence = mean(Cxy);
            coherence_matrix(i, j) = avg_coherence;
            coherence_matrix(j, i) = avg_coherence;  % Symmetric
        end
    end
    % Average coherence across all unique pairs
    spectral_coherence_index = mean(coherence_matrix(upper_triangle_indices));
    
    fprintf('  Global synchrony index: %.4f\n', global_synchrony_index);
    fprintf('  Spectral coherence index: %.4f\n', spectral_coherence_index);
    fprintf('  Significant connections: %d/%d (%.1f%%)\n', ...
            n_significant, numel(unique_correlations), percent_significant);
    fprintf('  Functional connectivity range: [%.3f, %.3f]\n', ...
            min(unique_correlations), max(unique_correlations));
    
    %% Save results
    fprintf('  Saving results...\n');
    
    % Load all original variables first
    original_data = load(fil);
    
    % Add new calculated variables
    original_data.burstiness = burstiness;
    original_data.firing_rate_hz = firing_rate_hz;
    original_data.P_spike = P_spike;
    original_data.R = R;
    original_data.R_pvalue = R_pvalue;
    original_data.R_significant = R_significant;
    original_data.coherence_matrix = coherence_matrix;
    original_data.global_synchrony_timeseries = global_synchrony_timeseries;
    original_data.global_synchrony_index = global_synchrony_index;
    original_data.spectral_coherence_index = spectral_coherence_index;
    original_data.params = params;
    original_data.analysis_metadata = struct(...
        'frame_rate', frame_rate, ...
        'max_lag', max_lag, ...
        'Nsamples', Nsamples, ...
        'n_shuffles', n_shuffles, ...
        'significance_alpha', significance_alpha, ...
        'n_rois', n_rois, ...
        'n_frames', n_frames, ...
        'total_time', total_time, ...
        'n_significant_connections', n_significant, ...
        'percent_significant_connections', percent_significant, ...
        'analysis_date', datetime('now'));
    
    % Save back to the same file (overwrite)
    save(fil, '-struct', 'original_data');
    
    fprintf('  File updated successfully.\n\n');
end

%% Summary
fprintf('========================================\n');
fprintf('Analysis complete!\n');
fprintf('Processed %d files.\n', numel(filelist));
fprintf('========================================\n');
fprintf('\nAdded variables to each file:\n');
fprintf('  - burstiness: Coefficient of variation of ISI for each ROI\n');
fprintf('  - firing_rate_hz: Firing rate in Hz for each ROI\n');
fprintf('  - P_spike: Spike probability matrix (ROIs x frames)\n');
fprintf('  - R: Functional connectivity matrix (ROIs x ROIs, signed cross-correlation)\n');
fprintf('  - R_pvalue: P-values for each connection (via spike train shuffling)\n');
fprintf('  - R_significant: Boolean mask for significant connections\n');
fprintf('  - coherence_matrix: Spectral coherence matrix (ROIs x ROIs)\n');
fprintf('  - global_synchrony_timeseries: Population synchrony at each time point\n');
fprintf('  - global_synchrony_index: Average global synchrony across time\n');
fprintf('  - spectral_coherence_index: Average spectral coherence across all pairs\n');
fprintf('  - params: Parameters used for analysis\n');
fprintf('  - analysis_metadata: Analysis settings and timestamp\n');
