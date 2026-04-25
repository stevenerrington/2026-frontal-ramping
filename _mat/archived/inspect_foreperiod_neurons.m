%% ========================================================================
% Foreperiod ramping: prevalence across cortical areas
% Steven Errington, 2026
% =========================================================================

linearity_analysis_table = readtable('linearity_analysis.csv');

%% ------------------------------------------------------------------------
% 1. Classify ramping neurons and append flag to analysis_log
% -------------------------------------------------------------------------

neuron_idx.pos_ramping = find(...
    linearity_analysis_table.flag  == 1  & ...
    linearity_analysis_table.slope >  0  & ...
    linearity_analysis_table.r2    >  0.8);

neuron_idx.neg_ramping = find(...
    linearity_analysis_table.flag  == 1  & ...
    linearity_analysis_table.slope <  0  & ...
    linearity_analysis_table.r2    >  0.8);

fprintf('Positive ramp neurons: %i\n', length(neuron_idx.pos_ramping))
fprintf('Negative ramp neurons: %i\n', length(neuron_idx.neg_ramping))

% Append ramping flag to log (1 = pos, -1 = neg, 0 = non-ramping)
analysis_log = spike_log;
analysis_log.ramping_fp_flag = zeros(size(analysis_log, 1), 1);
analysis_log.ramping_fp_flag(neuron_idx.pos_ramping) =  1;
analysis_log.ramping_fp_flag(neuron_idx.neg_ramping) = -1;

%% ------------------------------------------------------------------------
% 2. Compute prevalence per area
% -------------------------------------------------------------------------

areas  = unique(analysis_log.area);
nAreas = length(areas);

prevalence_table = table(...
    areas,            ...
    zeros(nAreas, 1), ...
    zeros(nAreas, 1), ...
    zeros(nAreas, 1), ...
    zeros(nAreas, 1), ...
    zeros(nAreas, 1), ...
    zeros(nAreas, 1), ...
    'VariableNames', {'area','n_total','n_pos','n_neg','n_nonramp','n_ramp','pct_ramp'});

for area_i = 1:nAreas
    area_idx  = strcmp(analysis_log.area, areas{area_i});
    flags     = analysis_log.ramping_fp_flag(area_idx);

    n_total   = sum(area_idx);
    n_pos     = sum(flags ==  1);
    n_neg     = sum(flags == -1);
    n_nonramp = sum(flags ==  0);
    n_ramp    = n_pos + n_neg;

    prevalence_table.n_total(area_i)   = n_total;
    prevalence_table.n_pos(area_i)     = n_pos;
    prevalence_table.n_neg(area_i)     = n_neg;
    prevalence_table.n_nonramp(area_i) = n_nonramp;
    prevalence_table.n_ramp(area_i)    = n_ramp;
    prevalence_table.pct_ramp(area_i)  = 100 * (n_ramp / n_total);
end

prevalence_table = sortrows(prevalence_table, 'pct_ramp', 'descend');

%% ------------------------------------------------------------------------
% 3. Reorder by cytoarchitectonic hierarchy and apply minimum n threshold
% -------------------------------------------------------------------------

area_order = {
    '24c',                                  ... % MFC
    '6DR',  '6DC',  '6VaVb',               ... % PMC
    '8B',   '8A',   '46d', '46df', '46v',  ... % dlPFC
    '44',   '45',                           ... % IFG
    '12r',  '12m',  '12o', '12l',          ... % vlPFC
    'AI',                                   ... % agranular insula
    '13l',  '13m',  '11ml',                ... % OFC
    'cd', 'pu',                   ... % striatum
    'AMG'};                                % amygdala

% Region colour map
region_colours = [
    0.6  0.2  0.6;   % MFC      - purple
    0.9  0.5  0.1;   % PMC      - orange
    0.2  0.6  0.9;   % dlPFC    - blue
    0.3  0.7  0.4;   % IFG      - green
    0.1  0.4  0.7;   % vlPFC    - dark blue
    0.7  0.7  0.3;   % AI       - yellow
    0.8  0.3  0.3;   % OFC      - red
    0.5  0.5  0.5;   % striatum - grey
    0.3  0.3  0.3];  % amygdala - dark grey

% Region index per area (in area_order sequence, before filtering)
region_map_full = [...
    ones(1,1)*1, ones(1,3)*2, ones(1,5)*3, ones(1,2)*4, ...
    ones(1,4)*5, ones(1,1)*6, ones(1,3)*7, ones(1,2)*8, ones(1,1)*9];

% Apply minimum n threshold then reorder
min_n              = 10;
pt                 = prevalence_table(prevalence_table.n_total >= min_n, :);
[~, order_idx]     = ismember(area_order, pt.area);
order_idx          = order_idx(order_idx > 0);
pt_ordered         = pt(order_idx, :);
region_map         = region_map_full(order_idx > 0);
nPlot              = height(pt_ordered);

%% ------------------------------------------------------------------------
% 4. Omnibus chi-square tests
% -------------------------------------------------------------------------

% --- 4a. Ramp vs non-ramp across all areas ---
cont_ramp       = [pt_ordered.n_ramp, pt_ordered.n_nonramp];
[chi2_omni_ramp, p_omni_ramp, df_omni_ramp] = chi2_contingency(cont_ramp);

fprintf('\n--- Omnibus: ramp vs non-ramp across areas ---\n')
fprintf('X2(%i) = %.2f, p = %.4f\n', df_omni_ramp, chi2_omni_ramp, p_omni_ramp)

% --- 4b. Pos vs neg ratio across areas (ramping neurons only) ---
ramp_rows       = pt_ordered.n_ramp > 0;
cont_posneg     = [pt_ordered.n_pos(ramp_rows), pt_ordered.n_neg(ramp_rows)];
[chi2_omni_pn, p_omni_pn, df_omni_pn] = chi2_contingency(cont_posneg);

fprintf('\n--- Omnibus: pos vs neg ratio across areas ---\n')
fprintf('X2(%i) = %.2f, p = %.4f\n', df_omni_pn, chi2_omni_pn, p_omni_pn)

%% ------------------------------------------------------------------------
% 5. Per-area OR (one-vs-rest): ramp vs non-ramp
% -------------------------------------------------------------------------

overall_ramp    = sum(pt_ordered.n_ramp);
overall_total   = sum(pt_ordered.n_total);
overall_nonramp = overall_total - overall_ramp;

chi2_stat  = nan(nPlot, 1);
chi2_p     = nan(nPlot, 1);
or_val     = nan(nPlot, 1);
or_ci_lo   = nan(nPlot, 1);
or_ci_hi   = nan(nPlot, 1);

for area_i = 1:nPlot

    % 2x2 contingency: [ramp, nonramp] x [this area, other areas]
    a = pt_ordered.n_ramp(area_i);
    b = pt_ordered.n_nonramp(area_i);
    c = overall_ramp    - a;
    d = overall_nonramp - b;

    % Chi-square with Yates correction
    contingency = [a, b; c, d];
    rowSums     = sum(contingency, 2);
    colSums     = sum(contingency, 1);
    N           = sum(contingency(:));
    expected    = (rowSums * colSums) / N;

    chi2_stat(area_i) = sum(sum((abs(contingency - expected) - 0.5).^2 ./ expected));
    chi2_p(area_i)    = 1 - chi2cdf(chi2_stat(area_i), 1);

    if any(expected(:) < 5)
        warning('Area %s: expected cell < 5, chi-square may be unreliable', ...
            pt_ordered.area{area_i})
    end

    % Odds ratio + 95% CI (Woolf; Haldane correction if zero cell)
    if any([a b c d] == 0)
        a = a + 0.5; b = b + 0.5; c = c + 0.5; d = d + 0.5;
    end
    log_or           = log((a * d) / (b * c));
    se_log_or        = sqrt(1/a + 1/b + 1/c + 1/d);
    or_val(area_i)   = exp(log_or);
    or_ci_lo(area_i) = exp(log_or - 1.96 * se_log_or);
    or_ci_hi(area_i) = exp(log_or + 1.96 * se_log_or);
end

% BH-FDR correction
chi2_p_fdr = bh_fdr(chi2_p);

% Append stats to ordered table
pt_ordered.chi2_stat  = chi2_stat;
pt_ordered.chi2_p_raw = chi2_p;
pt_ordered.chi2_p_fdr = chi2_p_fdr;
pt_ordered.sig_fdr    = chi2_p_fdr < 0.05;
pt_ordered.odds_ratio = or_val;
pt_ordered.or_ci_lo   = or_ci_lo;
pt_ordered.or_ci_hi   = or_ci_hi;

fprintf('\n--- Per-area OR results (ramp vs non-ramp, FDR corrected) ---\n')
disp(pt_ordered(:, {'area','n_total','n_ramp','pct_ramp', ...
    'odds_ratio','or_ci_lo','or_ci_hi','chi2_p_raw','chi2_p_fdr','sig_fdr'}))

%% ------------------------------------------------------------------------
% 6. Three-panel figure
% -------------------------------------------------------------------------

pos_pct = 100 * pt_ordered.n_pos ./ pt_ordered.n_total;
neg_pct = 100 * pt_ordered.n_neg ./ pt_ordered.n_total;

figure('Position', [100 100 900 850]);

% Omnibus p annotation strings
omni_ramp_str = sprintf('Omnibus ramp vs non-ramp: X2(%i)=%.1f, p=%.4f', ...
    df_omni_ramp, chi2_omni_ramp, p_omni_ramp);
omni_pn_str   = sprintf('Omnibus pos vs neg: X2(%i)=%.1f, p=%.4f', ...
    df_omni_pn, chi2_omni_pn, p_omni_pn);

% --- Panels 1 & 2 combined: stacked pos/neg ramp prevalence ---
subplot(3,1,[1 2]); hold on
title({omni_ramp_str, omni_pn_str, ...
    'Foreperiod ramping prevalence across frontal areas'}, 'FontSize', 8)

for area_i = 1:nPlot
    col = region_colours(region_map(area_i), :);
    
    % Positive ramp (upward from 0)
    bar(area_i, pos_pct(area_i), ...
        'FaceColor', col, ...
        'EdgeColor', 'none', 'BarWidth', 0.7)
    
    % Negative ramp (downward from 0)
    bar(area_i, -neg_pct(area_i), ...
        'FaceColor', col * 0.5, ...    % darker shade for negative
        'EdgeColor', 'none', 'BarWidth', 0.7)
    
end

% Zero line
yline(0, '-k', 'LineWidth', 0.8)

set(gca, 'XTick', 1:nPlot, 'XTickLabel', {}, 'FontSize', 9)
ylabel('Ramp prevalence (%)')
ylim([-30 30])
xlim([0 nPlot + 0.5])

% Manual y-axis labels to show positive values on both sides
yticks(-30:10:30)
yticklabels(arrayfun(@(x) num2str(abs(x)), -30:10:30, 'UniformOutput', false))

% Axis annotations
text(-0.8, 15, 'Positive', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'FontSize', 8)
text(-0.8, -15, 'Negative', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'FontSize', 8)

box off
add_region_dividers(region_map)

% --- Panel 3: Odds ratio forest plot (ramp vs non-ramp, per area) ---
subplot(3,1,3); hold on
yline(1, '--k', 'LineWidth', 0.8)

for area_i = 1:nPlot
    col   = region_colours(region_map(area_i), :);
    or    = pt_ordered.odds_ratio(area_i);
    ci_lo = pt_ordered.or_ci_lo(area_i);
    ci_hi = pt_ordered.or_ci_hi(area_i);

    % CI whisker
    plot([area_i area_i], [ci_lo ci_hi], '-', 'Color', col, 'LineWidth', 1.5)

    % Point: filled = FDR significant, open = not
    if pt_ordered.chi2_p_fdr(area_i) < 0.05
        scatter(area_i, or, 60, col, 'filled', 'MarkerEdgeColor', 'none')
    else
        scatter(area_i, or, 60, col, 'MarkerEdgeColor', col, 'LineWidth', 1.2)
    end

    % Significance markers above/below bars
    marker = get_sig_marker(pt_ordered.chi2_p_fdr(area_i));
    if ~isempty(marker)
        text(area_i, 6, marker, ...
            'HorizontalAlignment', 'center', 'FontSize', 9)
    end

    % n label
    text(area_i, 0.15, sprintf('n=%i', pt_ordered.n_total(area_i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 6, 'Rotation', 0)
end

set(gca, 'XTick', 1:nPlot, ...
         'XTickLabel', pt_ordered.area, ...
         'XTickLabelRotation', 45, ...
         'YScale', 'log', ...
         'FontSize', 9)
ylabel('Odds ratio vs rest (log scale)')
xlim([0 nPlot + 0.5]); ylim([0.05 20])
box off; add_region_dividers(region_map)

%% ========================================================================
% Helper functions
% =========================================================================

function plot_prevalence_bars(pct, pt_ordered, region_colours, region_map, nPlot)
% Coloured bars with FDR significance markers
    for area_i = 1:nPlot
        bar(area_i, pct(area_i), ...
            'FaceColor', region_colours(region_map(area_i), :), ...
            'EdgeColor', 'none', 'BarWidth', 0.7)
        marker = get_sig_marker(pt_ordered.chi2_p_fdr(area_i));
        if ~isempty(marker)
            text(area_i, pct(area_i) + 1, marker, ...
                'HorizontalAlignment', 'center', 'FontSize', 9)
        end
    end
end

function marker = get_sig_marker(p)
% Significance asterisk string
    if p < 0.001
        marker = '***';
    elseif p < 0.01
        marker = '**';
    elseif p < 0.05
        marker = '*';
    else
        marker = '';
    end
end

function add_region_dividers(region_map)
% Dashed vertical lines at region boundaries
    region_bounds = find(diff(region_map));
    for rb = region_bounds
        xline(rb + 0.5, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8)
    end
end

function [chi2_stat, p_val, df] = chi2_contingency(cont)
% Chi-square test on an arbitrary contingency matrix (no toolbox required)
    rowSums  = sum(cont, 2);
    colSums  = sum(cont, 1);
    N        = sum(cont(:));
    expected = (rowSums * colSums) / N;
    chi2_stat = sum(sum((cont - expected).^2 ./ expected));
    df        = (size(cont, 1) - 1) * (size(cont, 2) - 1);
    p_val     = 1 - chi2cdf(chi2_stat, df);
end

function p_fdr = bh_fdr(p)
% Benjamini-Hochberg FDR correction (no toolbox required)
    n = length(p);
    [p_sorted, sort_idx] = sort(p(:));
    p_fdr_sorted = p_sorted .* n ./ (1:n)';
    for k = n-1:-1:1
        p_fdr_sorted(k) = min(p_fdr_sorted(k), p_fdr_sorted(k+1));
    end
    p_fdr           = nan(size(p));
    p_fdr(sort_idx) = min(p_fdr_sorted, 1);
end