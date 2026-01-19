%% Example Visualization Script for Calcium Imaging Analysis
% This script demonstrates how to visualize the outputs from analyze_calcium_imaging.m
% It creates various plots including raster plots, connectivity matrices, and synchrony analysis

%% Clear workspace
clear;
clc;

%% Load analyzed data
fprintf('Loading analyzed data...\n');
fprintf('Please select the analyzed .mat file (e.g., STAR_Protoc_example_mcor.mat)...\n');
[filename, filepath] = uigetfile('*.mat', 'Select analyzed .mat file (e.g., STAR_Protoc_example_mcor.mat)');

if isequal(filename, 0)
    error('No file selected. Exiting.');
end

data_file = fullfile(filepath, filename);
data = load(data_file);
fprintf('Data loaded successfully!\n\n');

%% Extract key variables
P_spike = data.P_spike;
R = data.R;
R_significant = data.R_significant;
burstiness = data.burstiness;
firing_rate_hz = data.firing_rate_hz;
global_synchrony_timeseries = data.global_synchrony_timeseries;

[n_rois, n_frames] = size(P_spike);
frame_period = data.analysis_metadata.frame_period;
time_vector = (0:n_frames-1) * frame_period;

fprintf('Dataset info:\n');
fprintf('  Number of ROIs: %d\n', n_rois);
fprintf('  Number of frames: %d\n', n_frames);
fprintf('  Total time: %.2f seconds\n', data.analysis_metadata.total_time);
fprintf('  Frame period: %.4f seconds\n\n', frame_period);

%% Figure 1: Raster Plot of Spike Probabilities
figure('Name', 'Spike Raster Plot', 'Position', [100, 100, 1200, 600]);

% Raster plot (binarized at threshold 0.5)
% Note: P_spike represents the probability that a spike occurs in each time bin
% Threshold of 0.5 means showing only frames where there is >=50% chance a spike occurs
% You can adjust this threshold (e.g., 0.3, 0.7) based on desired stringency
subplot(2, 1, 1);
imagesc(time_vector, 1:n_rois, P_spike > 0.5);
colormap(gca, [1 1 1; 0 0 0]);  % White background, black spikes
xlabel('Time (seconds)');
ylabel('ROI #');
title('Spike Raster Plot (P_{spike} > 0.5)');
set(gca, 'YDir', 'normal');

% Continuous spike probability heatmap
subplot(2, 1, 2);
imagesc(time_vector, 1:n_rois, P_spike);
colormap(gca, hot);
colorbar;
xlabel('Time (seconds)');
ylabel('ROI #');
title('Spike Probability Heatmap');
set(gca, 'YDir', 'normal');

%% Figure 2: Population Activity and Synchrony
figure('Name', 'Population Activity', 'Position', [150, 150, 1200, 800]);

% Population firing rate over time
subplot(3, 1, 1);
population_rate = sum(P_spike, 1);
plot(time_vector, population_rate, 'k', 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Population Activity (spikes/frame)');
title('Population Firing Rate Over Time');
grid on;

% Global synchrony time series
subplot(3, 1, 2);
plot(time_vector, global_synchrony_timeseries, 'b', 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Synchrony Index');
title(sprintf('Global Population Synchrony (Mean = %.4f)', data.global_synchrony_index));
grid on;

% Cross-correlation of population activity with synchrony
subplot(3, 1, 3);
yyaxis left;
plot(time_vector, population_rate, 'k', 'LineWidth', 1);
ylabel('Population Activity');
yyaxis right;
plot(time_vector, global_synchrony_timeseries, 'b', 'LineWidth', 1);
ylabel('Synchrony Index');
xlabel('Time (seconds)');
title('Population Activity vs. Synchrony');
grid on;

%% Figure 3: Functional Connectivity Matrix
figure('Name', 'Functional Connectivity', 'Position', [200, 200, 1400, 500]);

% Raw connectivity matrix
subplot(1, 3, 1);
imagesc(R);
colormap(gca, jet);
colorbar;
caxis([-1, 1]);
xlabel('ROI #');
ylabel('ROI #');
title('Functional Connectivity Matrix (R)');
axis square;

% Significant connections only
subplot(1, 3, 2);
R_sig_only = R .* R_significant;
R_sig_only(~R_significant) = NaN;
imagesc(R_sig_only);
colormap(gca, jet);
colorbar;
caxis([-1, 1]);
xlabel('ROI #');
ylabel('ROI #');
title(sprintf('Significant Connections (p < %.2f)', data.analysis_metadata.significance_alpha));
axis square;

% Connection strength histogram
subplot(1, 3, 3);
upper_tri = triu(true(n_rois), 1);
all_connections = R(upper_tri);
sig_connections = R(upper_tri & R_significant);

histogram(all_connections, 30, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on;
histogram(sig_connections, 30, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('Correlation Coefficient');
ylabel('Count');
title('Distribution of Functional Connections');
legend({'All connections', 'Significant connections'}, 'Location', 'northwest');
grid on;

%% Figure 4: Single Cell Properties
figure('Name', 'Single Cell Properties', 'Position', [250, 250, 1200, 800]);

% Firing rate distribution
subplot(2, 2, 1);
histogram(firing_rate_hz, 20, 'FaceColor', 'b', 'EdgeColor', 'none');
xlabel('Firing Rate (Hz)');
ylabel('Count');
title(sprintf('Firing Rate Distribution (Mean = %.3f Hz)', mean(firing_rate_hz)));
grid on;

% Burstiness distribution
subplot(2, 2, 2);
histogram(burstiness, 20, 'FaceColor', 'g', 'EdgeColor', 'none');
xlabel('Burstiness (CV of ISI)');
ylabel('Count');
title(sprintf('Burstiness Distribution (Mean = %.3f)', mean(burstiness)));
grid on;

% Firing rate vs. Burstiness
subplot(2, 2, 3);
scatter(firing_rate_hz, burstiness, 50, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Firing Rate (Hz)');
ylabel('Burstiness (CV of ISI)');
title('Firing Rate vs. Burstiness');
grid on;

% Network connectivity per ROI
subplot(2, 2, 4);
connectivity_degree = sum(R_significant, 2) - 1;  % Subtract self-connection
bar(1:n_rois, connectivity_degree, 'FaceColor', [0.8 0.4 0.4]);
xlabel('ROI #');
ylabel('Number of Significant Connections');
title(sprintf('Network Degree (Mean = %.1f connections)', mean(connectivity_degree)));
grid on;

%% Figure 5: Example Single ROI Traces
figure('Name', 'Example ROI Traces', 'Position', [300, 300, 1200, 800]);

% Select up to 6 example ROIs (evenly spaced or most active)
[~, sorted_idx] = sort(firing_rate_hz, 'descend');
n_examples = min(6, n_rois);
example_rois = sorted_idx(1:n_examples);

for i = 1:n_examples
    subplot(n_examples, 1, i);
    roi_idx = example_rois(i);
    
    % Plot spike probability
    plot(time_vector, P_spike(roi_idx, :), 'k', 'LineWidth', 1);
    ylabel(sprintf('ROI %d', roi_idx));
    title(sprintf('ROI %d: FR = %.3f Hz, Burst = %.2f, Connections = %d', ...
        roi_idx, firing_rate_hz(roi_idx), burstiness(roi_idx), ...
        sum(R_significant(roi_idx, :)) - 1));
    grid on;
    
    if i == n_examples
        xlabel('Time (seconds)');
    else
        set(gca, 'XTickLabel', []);
    end
end

%% Print Summary Statistics
fprintf('========================================\n');
fprintf('ANALYSIS SUMMARY\n');
fprintf('========================================\n');
fprintf('Network Properties:\n');
fprintf('  Total ROIs: %d\n', n_rois);
fprintf('  Recording duration: %.2f seconds\n', data.analysis_metadata.total_time);
fprintf('  Frame period: %.4f seconds\n\n', frame_period);

fprintf('Firing Statistics:\n');
fprintf('  Mean firing rate: %.3f ± %.3f Hz\n', mean(firing_rate_hz), std(firing_rate_hz));
fprintf('  Firing rate range: [%.3f, %.3f] Hz\n', min(firing_rate_hz), max(firing_rate_hz));
fprintf('  Mean burstiness: %.3f ± %.3f\n', mean(burstiness), std(burstiness));
fprintf('  Burstiness range: [%.3f, %.3f]\n\n', min(burstiness), max(burstiness));

fprintf('Network Connectivity:\n');
fprintf('  Significant connections: %d/%d (%.1f%%)\n', ...
    data.analysis_metadata.n_significant_connections, ...
    numel(all_connections), ...
    data.analysis_metadata.percent_significant_connections);
fprintf('  Mean correlation: %.3f\n', mean(all_connections));
fprintf('  Mean significant correlation: %.3f\n', mean(sig_connections));
fprintf('  Mean network degree: %.1f connections/ROI\n\n', mean(connectivity_degree));

fprintf('Synchronization Metrics:\n');
fprintf('  Global synchrony index: %.4f\n', data.global_synchrony_index);
fprintf('  Spectral coherence index: %.4f\n', data.spectral_coherence_index);
fprintf('========================================\n');

fprintf('\nVisualization complete! All figures generated.\n');
