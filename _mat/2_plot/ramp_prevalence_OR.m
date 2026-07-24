%% ========================================================================
% Stacked bar + OR forest plot: foreperiod and reward ramping per area
% Top panel:    stacked bars — positive ramp up, negative ramp down
% Bottom panel: odds ratio forest plot, FDR corrected
% Self-contained — includes area ordering, region colours, filtering
% Steven Errington, 2026
% =========================================================================

%% ------------------------------------------------------------------------
% 1. Area order, region colours, and region map
% -------------------------------------------------------------------------

area_order = {
    '24c',                                  ... % MFC
    '6DR',  '6DC',  '6VaVb',               ... % PMC
    '8B',   '8A',   '46d', '46df', '46v',  ... % dlPFC
    '44',   '45',                           ... % IFG
    '12r',  '12m',  '12o', '12l',          ... % vlPFC
    'AI',                                   ... % agranular insula
    '13l',  '13m',  '11ml',                ... % OFC
    'cd',   'pu',                           ... % striatum
    'AMG'};                                 % amygdala

region_colours = [
    0.6  0.2  0.6;   % 1  MFC
    0.9  0.5  0.1;   % 2  PMC
    0.2  0.6  0.9;   % 3  dlPFC
    0.3  0.7  0.4;   % 4  IFG
    0.1  0.4  0.7;   % 5  vlPFC
    0.7  0.7  0.3;   % 6  AI
    0.8  0.3  0.3;   % 7  OFC
    0.5  0.5  0.5;   % 8  striatum
    0.3  0.3  0.3];  % 9  amygdala

region_map = [...
    1, ...
    2,  2,  2, ...
    3,  3,  3,  3,  3, ...
    4,  4, ...
    5,  5,  5,  5, ...
    6, ...
    7,  7,  7, ...
    8,  8, ...
    9];

%% ------------------------------------------------------------------------
% 2. Filter to areas present in analysis_log with >= MIN_N neurons
% -------------------------------------------------------------------------

MIN_N = 10;

analysis_log.area = strtrim(analysis_log.area);

[area_names_found, ~, area_ic] = unique(analysis_log.area);
area_n         = accumarray(area_ic, 1);
areas_with_n   = area_names_found(area_n >= MIN_N);

in_data         = ismember(area_order, areas_with_n);
area_order_plot = area_order(in_data);
region_map_plot = region_map(in_data);
nPlot           = numel(area_order_plot);

fprintf('Areas included: %i of %i\n', nPlot, numel(area_order));

%% ------------------------------------------------------------------------
% 3. Epoch definitions
% -------------------------------------------------------------------------

epochs = struct( ...
    'label', {'Foreperiod',      'Reward'}, ...
    'flag',  {'ramping_fp_flag', 'ramping_rwd_flag'});

nEpochs  = numel(epochs);
bar_width = 0.35;
x_centres = (1:nPlot);
x_offsets = [-bar_width/2, bar_width/2];
x_jitter  = [-0.18, 0.18];

%% ------------------------------------------------------------------------
% 4. COLOUR FIX ONLY (NEW)
% ------------------------------------------------------------------------

blend_to_white = @(c, f) c*f + (1 - f);

fp_scale  = 0.35;   % Foreperiod = lighter
rwd_scale = 0.85;   % Reward = darker

%% ------------------------------------------------------------------------
% 5. Compute per-area stats: prevalence + OR
% (UNCHANGED)
% -------------------------------------------------------------------------

pct_pos = nan(nPlot, nEpochs);
pct_neg = nan(nPlot, nEpochs);
or_val  = nan(nPlot, nEpochs);
ci_lo   = nan(nPlot, nEpochs);
ci_hi   = nan(nPlot, nEpochs);
chi2_p  = nan(nPlot, nEpochs);

for ep_i = 1:nEpochs

    flag_col        = epochs(ep_i).flag;
    flags           = analysis_log.(flag_col);
    overall_ramp    = sum(abs(flags) == 1);
    overall_total   = height(analysis_log);
    overall_nonramp = overall_total - overall_ramp;

    for area_i = 1:nPlot

        area_mask = strcmp(analysis_log.area, area_order_plot{area_i});
        n_total   = sum(area_mask);
        area_flags = flags(area_mask);

        n_pos = sum(area_flags ==  1);
        n_neg = sum(area_flags == -1);

        pct_pos(area_i, ep_i) =  100 * n_pos / n_total;
        pct_neg(area_i, ep_i) = -100 * n_neg / n_total;

        a = sum(abs(area_flags) == 1);
        b = sum(area_flags == 0);
        c = overall_ramp - a;
        d = overall_nonramp - b;

        cont = [a b; c d];

        expected = (sum(cont,2) * sum(cont,1)) / sum(cont(:));
        chi2 = sum((abs(cont - expected) - 0.5).^2 ./ expected, 'all');
        chi2_p(area_i, ep_i) = 1 - chi2cdf(chi2, 1);

        if any([a b c d] == 0)
            a = a+0.5; b = b+0.5; c = c+0.5; d = d+0.5;
        end

        log_or = log((a*d)/(b*c));
        se     = sqrt(1/a + 1/b + 1/c + 1/d);

        or_val(area_i, ep_i) = exp(log_or);
        ci_lo(area_i, ep_i)  = exp(log_or - 1.96*se);
        ci_hi(area_i, ep_i)  = exp(log_or + 1.96*se);
    end
end

% BH-FDR
p_flat = chi2_p(:);
valid  = ~isnan(p_flat);
p_fdr  = nan(size(p_flat));
p_fdr(valid) = bh_fdr(p_flat(valid));
p_fdr = reshape(p_fdr, nPlot, nEpochs);

%% ------------------------------------------------------------------------
% 6. FIGURE
% -------------------------------------------------------------------------

fig_w = max(1000, nPlot * 46);
figure('Position', [50 50 fig_w 720], 'Renderer', 'painters');

ax1 = subplot(3,1,[1 2]);
hold(ax1,'on');

reg_bounds = find(diff(region_map_plot));
for rb = reg_bounds
    xline(ax1, rb+0.5, '--', 'Color',[0.75 0.75 0.75]);
end

yline(ax1,0,'k');

%% ----------------- STACKED BARS (ONLY CHANGE IS COLOUR) -----------------
for ep_i = 1:nEpochs

    for area_i = 1:nPlot

        base_col = region_colours(region_map_plot(area_i), :);

        if ep_i == 1
            col = blend_to_white(base_col, fp_scale);
        else
            col = blend_to_white(base_col, rwd_scale);
        end

        x_pos = x_centres + x_offsets(ep_i);

        bar(ax1, x_pos(area_i), pct_pos(area_i, ep_i), bar_width, ...
            'FaceColor', col, 'EdgeColor', 'none');

        bar(ax1, x_pos(area_i), pct_neg(area_i, ep_i), bar_width, ...
            'FaceColor', col, 'EdgeColor', 'none');
    end
end

set(ax1, 'XTick', []);
ylabel(ax1, 'Ramp prevalence (%)');

ylim_val = 25;
ylim(ax1, [-ylim_val ylim_val]);

ytk = -ylim_val:5:ylim_val;
set(ax1,'YTick',ytk,'YTickLabel',arrayfun(@(v) num2str(abs(v)),ytk,'UniformOutput',false));

%% ----------------- FOREST PLOT (ONLY COLOUR CHANGED) -----------------
ax2 = subplot(3,1,3);
hold(ax2,'on');

yline(ax2,1,'k--');

for rb = reg_bounds
    xline(ax2, rb+0.5, '--', 'Color',[0.75 0.75 0.75]);
end

for ep_i = 1:nEpochs
    for area_i = 1:nPlot

        x = area_i + x_jitter(ep_i);

        base_col = region_colours(region_map_plot(area_i), :);

        if ep_i == 1
            col = blend_to_white(base_col, fp_scale);
        else
            col = blend_to_white(base_col, rwd_scale);
        end

        or = or_val(area_i,ep_i);
        lo = ci_lo(area_i,ep_i);
        hi = ci_hi(area_i,ep_i);
        p  = p_fdr(area_i,ep_i);

        if isnan(or), continue; end

        plot([x x],[lo hi],'Color',col,'LineWidth',1.4);

        if p < 0.05
            scatter(x,or,50,col,'filled','MarkerEdgeColor','none');
        else
            scatter(x,or,50,col,'MarkerEdgeColor',col,'LineWidth',1.2);
        end

        if p < 0.05
            text(x, 5, sig_marker_local(p), ...
                'HorizontalAlignment','center','FontSize',8);
        end
    end
end

set(ax2,'YScale','log',...
    'XTick',1:nPlot,...
    'XTickLabel',area_order_plot,...
    'XTickLabelRotation',90,...
    'YLim', [0 10]);

ylabel(ax2,'Odds ratio');

%% ------------------------------------------------------------------------
% helpers (UNCHANGED)
% -------------------------------------------------------------------------

function m = sig_marker_local(p)
if isnan(p) || p >= 0.05, m = ''; return; end
if p < 0.001, m = '***'; return; end
if p < 0.01, m = '**'; return; end
m = '*';
end

