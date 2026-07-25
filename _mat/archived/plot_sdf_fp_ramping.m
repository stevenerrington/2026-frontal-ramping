%% ========================================================================
% Foreperiod ramping: average SDF per area, split by foreperiod condition
% Separate figures for positive and negative ramp neurons
% Steven Errington, 2026
% =========================================================================

%% ------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------

time_vec   = -2000:2000;
plot_win   = time_vec >= -1500 & time_vec <= 250;
time_plot  = time_vec(plot_win);

cond_labels = {'All', 'Short', 'Medium', 'Long'};
cond_colors = [0.5 0.5 0.5;   % all    - grey
               0.2 0.6 0.9;   % short  - light blue
               0.1 0.4 0.7;   % medium - mid blue
               0.0 0.2 0.5];  % long   - dark blue

fp_times  = [-1300, -1000, -700];
areas     = unique(analysis_log.area);
nAreas    = length(areas);
nCols     = 4;
nRows     = ceil(nAreas / nCols);

% Ramp populations to plot
ramp_types = struct(...
    'label', {'Positive ramp', 'Negative ramp'}, ...
    'flag',  {1, -1});

%% ------------------------------------------------------------------------
% One figure per ramp type
% -------------------------------------------------------------------------

for ramp_i = 1:length(ramp_types)

    flag  = ramp_types(ramp_i).flag;
    label = ramp_types(ramp_i).label;

    figure('Renderer', 'painters', 'Position', [50 50 700 nRows * 250]);
    sgtitle(sprintf('Foreperiod SDF — %s neurons', label), 'FontSize', 13)

    for area_i = 1:nAreas

        % Neurons in this area with the relevant ramp flag
        area_idx = strcmp(analysis_log.area, areas{area_i}) & ...
                   analysis_log.ramping_fp_flag == flag;
        n_neurons = sum(area_idx);

        subplot(nRows, nCols, area_i); hold on

        % Title regardless of whether neurons exist
        title(sprintf('%s  (n=%i)', areas{area_i}, n_neurons), 'FontSize', 10)
        xlabel('Time from stim on (ms)', 'FontSize', 8)
        ylabel('Firing rate (z)', 'FontSize', 8)
        xlim([-1500 250])
        box off

        % Skip plotting if no neurons in this area
        if n_neurons == 0
            text(mean([-1500 250]), 0, 'No neurons', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'Color', [0.6 0.6 0.6])
            continue
        end

        area_data = fixation_array(area_idx, plot_win, :);

        for cond_i = 2:4
            cond_data = squeeze(area_data(:, :, cond_i));  % n_neurons x time

            % Handle single neuron edge case
            if n_neurons == 1
                mu  = cond_data;
                sem = zeros(size(mu));
            else
                mu  = nanmean(cond_data, 1);
                sem = nanstd(cond_data, 0, 1) / sqrt(n_neurons);
            end

            % SEM shading
            % fill([time_plot, fliplr(time_plot)], ...
            %      [mu + sem, fliplr(mu - sem)], ...
            %      cond_colors(cond_i, :), ...
            %      'FaceAlpha', 0.15, 'EdgeColor', 'none')

            % Mean line
            plot(time_plot, mu, ...
                'Color', cond_colors(cond_i, :), ...
                'LineWidth', 1.5, ...
                'DisplayName', cond_labels{cond_i})
        end

        % Reference lines
        for fp = fp_times
            xline(fp, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8)
        end
        xline(0, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.0)

        % Legend on first subplot only
        if area_i == 1
            legend(cond_labels, 'Location', 'northwest', 'FontSize', 7)
        end

    end
end