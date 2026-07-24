%% ========================================================================
% Venn diagram: foreperiod vs reward ramping overlap
% Collapsed across all areas, positive and negative combined
% =========================================================================

% ── Classify neurons ──────────────────────────────────────────────────────
% Any ramping = abs(flag) == 1, regardless of direction
ramp_fp  = abs(analysis_log.ramping_fp_flag)  == 1;
ramp_rwd = abs(analysis_log.ramping_rwd_flag) == 1;

n_total    = height(analysis_log);
n_fp_only  = sum( ramp_fp & ~ramp_rwd);
n_rwd_only = sum(~ramp_fp &  ramp_rwd);
n_both     = sum( ramp_fp &  ramp_rwd);
n_neither  = sum(~ramp_fp & ~ramp_rwd);

fprintf('\n--- Foreperiod vs Reward ramping overlap ---\n')
fprintf('Total neurons      : %i\n',   n_total)
fprintf('Foreperiod only    : %i (%.1f%%)\n', n_fp_only,  100*n_fp_only/n_total)
fprintf('Reward only        : %i (%.1f%%)\n', n_rwd_only, 100*n_rwd_only/n_total)
fprintf('Both               : %i (%.1f%%)\n', n_both,     100*n_both/n_total)
fprintf('Neither            : %i (%.1f%%)\n', n_neither,  100*n_neither/n_total)

% ── Fisher exact test ─────────────────────────────────────────────────────
% 2x2: [ramps_fp ramps_rwd] x [yes no]
cont = [n_both,     n_fp_only; ...
        n_rwd_only, n_neither];
[~, p_fisher] = fishertest(cont);
fprintf('Fisher exact test  : p = %.4f\n', p_fisher)

% ── Draw Venn ─────────────────────────────────────────────────────────────
% Circle geometry: two overlapping circles, overlap area proportional to n_both
% We set circle separation to give visually clear but honest overlap

theta   = linspace(0, 2*pi, 300);
r       = 1.0;           % circle radius (arbitrary units)
sep     = 1.1;           % centre-to-centre distance (< 2r = overlap)

cx_fp   = -sep/2;        % centre x for foreperiod circle
cx_rwd  =  sep/2;        % centre x for reward circle
cy      = 0;

% Circle coordinates
xc_fp  = cx_fp  + r * cos(theta);
yc_fp  = cy     + r * sin(theta);
xc_rwd = cx_rwd + r * cos(theta);
yc_rwd = cy     + r * sin(theta);

% Colours
col_fp  = [0.20 0.60 0.90];   % blue  — foreperiod
col_rwd = [0.90 0.40 0.20];   % red   — reward
col_both = (col_fp + col_rwd) / 2;  % blended for overlap

figure('Position', [100 100 520 420], 'Renderer', 'painters');
ax = axes; hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');

% Draw filled circles (foreperiod behind reward)
patch(ax, xc_fp,  yc_fp,  col_fp,  'FaceAlpha', 0.45, 'EdgeColor', col_fp,  'LineWidth', 1.5);
patch(ax, xc_rwd, yc_rwd, col_rwd, 'FaceAlpha', 0.45, 'EdgeColor', col_rwd, 'LineWidth', 1.5);

% ── Count labels ─────────────────────────────────────────────────────────
pct = @(n) sprintf('%i\n(%.1f%%)', n, 100*n/n_total);

% FP only: left of overlap
text(ax, cx_fp - r*0.42, cy, pct(n_fp_only), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, ...
    'FontWeight', 'bold', 'Color', col_fp * 0.6);

% Reward only: right of overlap
text(ax, cx_rwd + r*0.42, cy, pct(n_rwd_only), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, ...
    'FontWeight', 'bold', 'Color', col_rwd * 0.6);

% Both: centre overlap
text(ax, 0, cy, pct(n_both), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, ...
    'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

% ── Circle labels ─────────────────────────────────────────────────────────
text(ax, cx_fp  - r*0.3, r*1.15, 'Foreperiod', ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', col_fp * 0.7);
text(ax, cx_rwd + r*0.3, r*1.15, 'Reward', ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', col_rwd * 0.7);

% ── Title with stats ──────────────────────────────────────────────────────
title(ax, sprintf('n = %i neurons   Fisher p = %.4f', n_total, p_fisher), ...
    'FontSize', 9);

% ── Neither label (outside both circles) ─────────────────────────────────
text(ax, 0, -r*1.4, sprintf('Neither: %i (%.1f%%)', n_neither, 100*n_neither/n_total), ...
    'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);

xlim(ax, [-r*2.2  r*2.2]);
ylim(ax, [-r*1.8  r*1.8]);

exportgraphics(gcf, fullfile(dirs.figures, 'venn_fp_rwd.pdf'), 'ContentType', 'vector');
fprintf('Saved: venn_fp_rwd.pdf\n')