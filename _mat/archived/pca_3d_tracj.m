%% ========================================================================
% pca_3d_trajectories.m
% 3D population trajectory plots in PC1 x PC2 x PC3 space
% -------------------------------------------------------------------------
% Run after jazayeri_pca_trajectories.m — requires traj{} and traj_norm{}
% in the workspace.
%
% Produces three figures:
%   1. Stim-aligned 3D trajectories (tests convergence at stim onset)
%   2. Time-normalised 3D trajectories (tests temporal scaling)
%   3. 2x2 panel: 3D + three 2D projections side by side
%
% Steven Errington, 2026
% =========================================================================

%% ------------------------------------------------------------------------
% Helper: draw a smooth arrow in 3D to show trajectory direction
% -------------------------------------------------------------------------

function draw_arrow3(ax, x, y, z, col, sz)
    % Draws an arrowhead at the end of a 3D trajectory segment
    % using a scatter marker — MATLAB's quiver3 is unreliable for this
    scatter3(ax, x, y, z, sz, col, 'filled', ...
        'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
end

%% ------------------------------------------------------------------------
% 0. Check prerequisites
% -------------------------------------------------------------------------

if ~exist('traj', 'var') || ~exist('traj_norm', 'var')
    error('Run jazayeri_pca_trajectories.m first to generate traj{} and traj_norm{}.');
end

nConds       = numel(traj);
fp_durations = [725, 1000, 1300];
cond_labels  = {'Short (725 ms)', 'Medium (1000 ms)', 'Long (1300 ms)'};
fp_colors    = [0.17 0.47 0.84;
                0.11 0.69 0.48;
                0.93 0.63 0.00];
N_NORM_PTS   = size(traj_norm{1}, 1);
t_norm_axis  = linspace(0, 1, N_NORM_PTS);

%% ========================================================================
% FIGURE 1: Stim-aligned 3D trajectories
% =========================================================================

fig1 = figure('Position', [50 50 620 540], 'Renderer', 'painters');
ax = axes(fig1);
hold(ax, 'on');
grid(ax, 'on');
ax.GridAlpha     = 0.15;
ax.GridLineStyle = ':';
view(ax, [-38 28]);   % initial viewing angle — user can rotate interactively

for c_i = 1:nConds
    tr  = traj{c_i};
    col = fp_colors(c_i, :);
    dur = fp_durations(c_i);

    % Find index of foreperiod start and stim onset
    fp_start_idx = find(t_full >= -dur, 1);
    t0_idx       = find(t_full >= 0, 1);

    % Full trajectory line
    plot3(ax, tr(:,1), tr(:,2), tr(:,3), '-', ...
        'Color', [col, 0.85], 'LineWidth', 2.5, ...
        'DisplayName', cond_labels{c_i});

    % Foreperiod start marker (open circle)
    scatter3(ax, tr(fp_start_idx,1), tr(fp_start_idx,2), tr(fp_start_idx,3), ...
        60, col, 'o', 'LineWidth', 1.8, 'HandleVisibility', 'off');

    % Stim onset marker (filled circle, black edge)
    scatter3(ax, tr(t0_idx,1), tr(t0_idx,2), tr(t0_idx,3), ...
        100, col, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    % Direction arrow: small dot slightly before end to show direction
    n_pts = size(tr, 1);
    mid   = round(fp_start_idx + (t0_idx - fp_start_idx) * 0.6);
    draw_arrow3(ax, tr(mid,1), tr(mid,2), tr(mid,3), col, 30);
end

% Stim onset convergence cloud (highlight where all conditions arrive at t=0)
t0_pts = cell2mat(cellfun(@(t) t(find(t_full>=0,1), :), traj, ...
    'UniformOutput', false));
scatter3(ax, t0_pts(:,1), t0_pts(:,2), t0_pts(:,3), 150, ...
    [0.85 0.15 0.15], 'd', 'filled', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1, 'DisplayName', 'Stim onset');

xlabel(ax, 'PC 1', 'FontSize', 10);
ylabel(ax, 'PC 2', 'FontSize', 10);
zlabel(ax, 'PC 3', 'FontSize', 10);
title(ax, 'Population trajectories — stim onset aligned', 'FontSize', 11);
legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 9);
set(ax, 'FontSize', 9, 'TickDir', 'out', 'Box', 'off');

% Drop lines to floor (aids 3D depth perception)
zlims = zlim(ax);
for c_i = 1:nConds
    t0_idx = find(t_full >= 0, 1);
    tr     = traj{c_i};
    plot3(ax, [tr(t0_idx,1) tr(t0_idx,1)], ...
              [tr(t0_idx,2) tr(t0_idx,2)], ...
              [tr(t0_idx,3) zlims(1)], ':', ...
        'Color', [fp_colors(c_i,:), 0.3], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end

exportgraphics(fig1, fullfile(dirs.figures, 'pca_3d_stim_aligned.pdf'), ...
    'ContentType', 'vector');
fprintf('Saved: pca_3d_stim_aligned.pdf\n');

%% ========================================================================
% FIGURE 2: Time-normalised 3D trajectories
% =========================================================================

fig2 = figure('Position', [100 100 620 540], 'Renderer', 'painters');
ax2 = axes(fig2);
hold(ax2, 'on');
grid(ax2, 'on');
ax2.GridAlpha     = 0.15;
ax2.GridLineStyle = ':';
view(ax2, [-38 28]);

for c_i = 1:nConds
    tr  = traj_norm{c_i};
    col = fp_colors(c_i, :);

    % Colour gradient along trajectory (light to full colour)
    n_pts = size(tr, 1);
    for seg_i = 1:n_pts-1
        alpha = 0.3 + 0.7 * (seg_i / n_pts);
        plot3(ax2, tr(seg_i:seg_i+1, 1), tr(seg_i:seg_i+1, 2), ...
              tr(seg_i:seg_i+1, 3), '-', ...
            'Color', [col * alpha + (1-alpha), 0.9], 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    end

    % Legend proxy
    plot3(ax2, nan, nan, nan, '-', 'Color', col, 'LineWidth', 2.5, ...
        'DisplayName', cond_labels{c_i});

    % Start (fixation onset, t_norm=0)
    scatter3(ax2, tr(1,1), tr(1,2), tr(1,3), 60, col, 'o', ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');

    % End (stim onset, t_norm=1)
    scatter3(ax2, tr(end,1), tr(end,2), tr(end,3), 100, col, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

xlabel(ax2, 'PC 1', 'FontSize', 10);
ylabel(ax2, 'PC 2', 'FontSize', 10);
zlabel(ax2, 'PC 3', 'FontSize', 10);
title(ax2, 'Time-normalised trajectories — temporal scaling test', 'FontSize', 11);
legend(ax2, 'Location', 'best', 'Box', 'off', 'FontSize', 9);
set(ax2, 'FontSize', 9, 'TickDir', 'out', 'Box', 'off');

% Annotation
text(ax2, 0.05, 0.95, ...
    'Overlap = temporal scaling supported', ...
    'Units', 'normalized', 'FontSize', 8, 'Color', [0.5 0.5 0.5], ...
    'VerticalAlignment', 'top');

exportgraphics(fig2, fullfile(dirs.figures, 'pca_3d_normalised.pdf'), ...
    'ContentType', 'vector');
fprintf('Saved: pca_3d_normalised.pdf\n');

%% ========================================================================
% FIGURE 3: Combined panel — 3D + three 2D projections
% =========================================================================

fig3 = figure('Position', [150 50 1000 450], 'Renderer', 'painters');

% ── Main 3D panel ────────────────────────────────────────────────────────
ax_3d = subplot(1, 4, 1);
hold(ax_3d, 'on');
grid(ax_3d, 'on');
ax_3d.GridAlpha = 0.15;
view(ax_3d, [-38 28]);

for c_i = 1:nConds
    tr  = traj{c_i};
    col = fp_colors(c_i,:);
    fp_start_idx = find(t_full >= -fp_durations(c_i), 1);
    t0_idx       = find(t_full >= 0, 1);

    plot3(ax_3d, tr(:,1), tr(:,2), tr(:,3), '-', ...
        'Color', [col 0.85], 'LineWidth', 2, ...
        'DisplayName', cond_labels{c_i});
    scatter3(ax_3d, tr(fp_start_idx,1), tr(fp_start_idx,2), ...
        tr(fp_start_idx,3), 40, col, 'o', 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    scatter3(ax_3d, tr(t0_idx,1), tr(t0_idx,2), tr(t0_idx,3), ...
        70, col, 'filled', 'MarkerEdgeColor','k','LineWidth',1, ...
        'HandleVisibility','off');
end

xlabel(ax_3d,'PC1','FontSize',8); ylabel(ax_3d,'PC2','FontSize',8);
zlabel(ax_3d,'PC3','FontSize',8);
title(ax_3d,'3D','FontSize',9);
legend(ax_3d,'Location','best','Box','off','FontSize',6);
set(ax_3d,'FontSize',7,'TickDir','out');

% ── PC1 vs PC2 ───────────────────────────────────────────────────────────
ax_12 = subplot(1, 4, 2);
hold(ax_12, 'on');
for c_i = 1:nConds
    tr  = traj{c_i};
    col = fp_colors(c_i,:);
    fp_start_idx = find(t_full >= -fp_durations(c_i), 1);
    t0_idx       = find(t_full >= 0, 1);
    plot(ax_12, tr(:,1), tr(:,2), '-', 'Color', [col 0.85], 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    scatter(ax_12, tr(fp_start_idx,1), tr(fp_start_idx,2), 40, col, 'o', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    scatter(ax_12, tr(t0_idx,1), tr(t0_idx,2), 70, col, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end
xlabel(ax_12,'PC1','FontSize',8); ylabel(ax_12,'PC2','FontSize',8);
title(ax_12,'PC1 vs PC2','FontSize',9);
box off; set(ax_12,'TickDir','out','FontSize',7);

% ── PC1 vs PC3 ───────────────────────────────────────────────────────────
ax_13 = subplot(1, 4, 3);
hold(ax_13, 'on');
for c_i = 1:nConds
    tr  = traj{c_i};
    col = fp_colors(c_i,:);
    fp_start_idx = find(t_full >= -fp_durations(c_i), 1);
    t0_idx       = find(t_full >= 0, 1);
    plot(ax_13, tr(:,1), tr(:,3), '-', 'Color', [col 0.85], 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    scatter(ax_13, tr(fp_start_idx,1), tr(fp_start_idx,3), 40, col, 'o', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    scatter(ax_13, tr(t0_idx,1), tr(t0_idx,3), 70, col, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end
xlabel(ax_13,'PC1','FontSize',8); ylabel(ax_13,'PC3','FontSize',8);
title(ax_13,'PC1 vs PC3','FontSize',9);
box off; set(ax_13,'TickDir','out','FontSize',7);

% ── PC2 vs PC3 ───────────────────────────────────────────────────────────
ax_23 = subplot(1, 4, 4);
hold(ax_23, 'on');
for c_i = 1:nConds
    tr  = traj{c_i};
    col = fp_colors(c_i,:);
    fp_start_idx = find(t_full >= -fp_durations(c_i), 1);
    t0_idx       = find(t_full >= 0, 1);
    plot(ax_23, tr(:,2), tr(:,3), '-', 'Color', [col 0.85], 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    scatter(ax_23, tr(fp_start_idx,2), tr(fp_start_idx,3), 40, col, 'o', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    scatter(ax_23, tr(t0_idx,2), tr(t0_idx,3), 70, col, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end
xlabel(ax_23,'PC2','FontSize',8); ylabel(ax_23,'PC3','FontSize',8);
title(ax_23,'PC2 vs PC3','FontSize',9);
box off; set(ax_23,'TickDir','out','FontSize',7);

% Shared legend annotation
annotation(fig3,'textbox',[0.01 0.02 0.98 0.06], ...
    'String', 'Open circle = fixation onset   Filled circle = stim onset', ...
    'EdgeColor','none','FontSize',8,'Color',[0.5 0.5 0.5], ...
    'HorizontalAlignment','center');

sgtitle('Population trajectories — PC space (stim aligned)', 'FontSize', 11);

exportgraphics(fig3, fullfile(dirs.figures, 'pca_3d_combined.pdf'), ...
    'ContentType', 'vector');
fprintf('Saved: pca_3d_combined.pdf\n');