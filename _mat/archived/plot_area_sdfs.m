function plot_area_sdfs(z_array, time_vec, analysis_log, area_order, ...
                        region_colours, varargin)
% PLOT_AREA_SDFS  Population-average SDF per area, faceted by cytoarchitectonic
%                 hierarchy, for positive and negative ramping neurons separately.
%
% Syntax
% ------
%   plot_area_sdfs(z_array, time_vec, analysis_log, area_order, region_colours)
%   plot_area_sdfs(__, Name, Value)
%
% Required inputs
% ---------------
%   z_array        : double, (nNeurons x nTimepoints x nConditions)
%                    Z-scored SDF array.
%   time_vec       : double, (1 x nTimepoints)
%                    Time vector in ms relative to alignment event.
%   analysis_log   : table
%                    Must contain columns: 'area', 'ramping_fp_flag'.
%                    ramping_fp_flag: 1 = positive ramp, -1 = negative, 0 = non.
%   area_order     : cell array of char
%                    Area labels in cytoarchitectonic order.
%   region_colours : double, (nRegions x 3)
%                    RGB colours for each region (matched to area_order via
%                    region_map — see below).
%
% Name-Value options
% ------------------
%   'ConditionLabels'  : cell array of char (default: {'All','Short','Med','Long'})
%   'ConditionColors'  : (nConds x 3) RGB — default: blues light to dark
%   'ConditionIdx'     : indices into z_array dim 3 to plot (default: all)
%   'XLim'            : [xmin xmax] in ms   (default: [time_vec(1) time_vec(end)])
%   'YLim'            : [ymin ymax] or []   (default: auto per subplot)
%   'VLines'          : ms values for vertical reference lines (default: [0])
%   'VLineStyle'      : line spec string    (default: 'k-')
%   'PlotSEM'         : logical, shade SEM  (default: true)
%   'MinN'            : minimum neurons to plot a panel (default: 5)
%   'NCols'           : number of subplot columns (default: 4)
%   'FigureTitle'     : char prefix for figure title (default: '')
%   'SavePath'        : char/string, full path to save PDF (default: '' = no save)
%   'RegionMap'       : int vector length = numel(area_order), maps each area
%                       to a row in region_colours. (default: ones, all same colour)
%
% Output
% ------
%   Two figures are created (positive ramp, negative ramp).
%
% Example
% -------
%   plot_area_sdfs(z_fixation_array, -2000:2000, analysis_log, ...
%       area_order, region_colours, ...
%       'ConditionLabels', {'All','Short','Medium','Long'}, ...
%       'ConditionIdx',    [1 2 3 4], ...
%       'XLim',            [-1500 250], ...
%       'VLines',          [-1300 -1000 -700 0], ...
%       'FigureTitle',     'Foreperiod SDF', ...
%       'SavePath',        fullfile(dirs.figures, 'foreperiod_sdf'));
%
% Steven Errington, 2026

% ── Parse inputs ─────────────────────────────────────────────────────────────
p = inputParser;
addRequired(p, 'z_array');
addRequired(p, 'time_vec');
addRequired(p, 'analysis_log');
addRequired(p, 'area_order');
addRequired(p, 'region_colours');

nConds_default = size(z_array, 3);
default_cond_labels = arrayfun(@(x) sprintf('Cond %i', x), ...
    1:nConds_default, 'UniformOutput', false);

% Default condition colours: blues light to dark
default_cond_colors = [
    0.36 0.76 0.94;   % light blue
    0.10 0.40 0.70;   % mid blue
    0.00 0.17 0.43;   % dark blue
    0.60 0.60 0.60;   % grey (fallback for extra conditions)
];
if nConds_default > size(default_cond_colors, 1)
    default_cond_colors = repmat(default_cond_colors(end,:), nConds_default, 1);
end

addParameter(p, 'ConditionLabels', default_cond_labels);
addParameter(p, 'ConditionColors', default_cond_colors);
addParameter(p, 'ConditionIdx',    1:nConds_default);
addParameter(p, 'XLim',           [time_vec(1) time_vec(end)]);
addParameter(p, 'YLim',           []);
addParameter(p, 'VLines',         0);
addParameter(p, 'VLineStyle',     'k-');
addParameter(p, 'PlotSEM',        true);
addParameter(p, 'MinN',           5);
addParameter(p, 'NCols',          4);
addParameter(p, 'FigureTitle',    '');
addParameter(p, 'SavePath',       '');
addParameter(p, 'RegionMap',      ones(1, numel(area_order)));

parse(p, z_array, time_vec, analysis_log, area_order, region_colours, varargin{:});
opts = p.Results;

cond_idx    = opts.ConditionIdx;
cond_labels = opts.ConditionLabels(cond_idx);
cond_colors = opts.ConditionColors(cond_idx, :);
nConds_plot = numel(cond_idx);
nAreas      = numel(area_order);
nCols       = opts.NCols;
nRows       = ceil(nAreas / nCols);

% Mask time vector to XLim for plotting
t_mask  = time_vec >= opts.XLim(1) & time_vec <= opts.XLim(2);
t_plot  = time_vec(t_mask);

% ── Strip whitespace from area column (defensive) ────────────────────────────
analysis_log.area = strtrim(analysis_log.area);

% ── Loop over ramp types ──────────────────────────────────────────────────────
ramp_types  = {  1,          -1         };
ramp_labels = {'Positive ramp', 'Negative ramp'};

for ramp_i = 1:2

    flag       = ramp_types{ramp_i};
    ramp_label = ramp_labels{ramp_i};

    fig = figure('Renderer', 'painters', ...
                 'Position', [50 50 nCols*280 nRows*220]);

    if ~isempty(opts.FigureTitle)
        sgtitle(sprintf('%s — %s neurons', opts.FigureTitle, ramp_label), ...
                'FontSize', 13);
    else
        sgtitle(sprintf('%s neurons', ramp_label), 'FontSize', 13);
    end

    for area_i = 1:nAreas

        area_name = area_order{area_i};
        col       = region_colours(opts.RegionMap(area_i), :);

        % Neuron indices: this area + this ramp type
        area_mask = strcmp(analysis_log.area, area_name) & ...
                    analysis_log.ramping_decision_flag == flag;
        n_neurons = sum(area_mask);

        ax = subplot(nRows, nCols, area_i);
        hold(ax, 'on');

        title(ax, sprintf('%s  (n=%i)', area_name, n_neurons), 'FontSize', 9);
        xlabel(ax, 'Time (ms)', 'FontSize', 7);
        ylabel(ax, 'FR (z)', 'FontSize', 7);
        set(ax, 'FontSize', 7, 'Box', 'off', 'TickDir', 'out');
        xlim(ax, opts.XLim);

        % Reference lines
        for vl = opts.VLines
            xline(ax, vl, opts.VLineStyle, 'LineWidth', 0.8, ...
                  'HandleVisibility', 'off');
        end

        % Skip if too few neurons
        if n_neurons < opts.MinN
            text(ax, mean(opts.XLim), 0, ...
                 sprintf('n < %i', opts.MinN), ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 8, 'Color', [0.7 0.7 0.7]);
            if ~isempty(opts.YLim)
                ylim(ax, opts.YLim);
            end
            continue
        end

        % Extract data for this area: (n_neurons x nTimepointsFull x nConds)
        area_data = z_array(area_mask, :, :);

        for cond_i = 1:nConds_plot

            ci     = cond_idx(cond_i);
            cdata  = squeeze(area_data(:, t_mask, ci));   % n_neurons x t_plot

            if n_neurons == 1
                mu  = cdata;
                sem = zeros(size(mu));
            else
                mu  = nanmean(cdata, 1);
                sem = nanstd(cdata, 0, 1) / sqrt(n_neurons);
            end

            c = cond_colors(cond_i, :);

            % SEM shading
            if opts.PlotSEM
                t_patch  = [t_plot, fliplr(t_plot)];
                y_patch  = [mu + sem, fliplr(mu - sem)];
                patch(ax, t_patch, y_patch, c, ...
                      'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
                      'HandleVisibility', 'off');
            end

            % Mean line
            plot(ax, t_plot, mu, '-', ...
                 'Color', c, 'LineWidth', 1.5, ...
                 'DisplayName', cond_labels{cond_i});
        end

        % Legend on first subplot only
        if area_i == 1
            legend(ax, 'Location', 'northwest', 'FontSize', 7, 'Box', 'off');
        end

        if ~isempty(opts.YLim)
            ylim(ax, opts.YLim);
        end

    end % area loop

    % Save if requested (append ramp type to filename)
    if ~isempty(opts.SavePath)
        ramp_suffix = strrep(lower(ramp_label), ' ', '_');
        save_file   = sprintf('%s_%s.pdf', opts.SavePath, ramp_suffix);
        exportgraphics(fig, save_file, 'ContentType', 'vector');
        fprintf('Saved: %s\n', save_file);
    end

end % ramp type loop

end