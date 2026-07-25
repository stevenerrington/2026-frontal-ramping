%% ========================================================================
% jazayeri_pca_trajectories.m
% Population trajectory analysis for elapsed time vs time-to-event coding
% -------------------------------------------------------------------------
% Implements the Jazayeri & Shadlen (2010) / Remington et al. approach:
%
%   1. Define PCA subspace on common window (-725 to 0 ms, all conditions)
%   2. Project full-length condition trajectories into that subspace
%   3. Plot trajectories aligned to stim onset (tests convergence)
%   4. Plot trajectories aligned to fixation onset (tests common origin)
%   5. Time-normalise and replot (tests temporal scaling)
%   6. Compute principal angles between condition subspaces (formal test)
%
% Required workspace variables:
%   z_fixation_array : nNeurons x nTime x nConds (z-scored SDFs)
%                      time_vec = -2000:2000, conds 2-4 = short/med/long
%   analysis_log     : table with 'area', 'ramping_fp_flag'
%
% Steven Errington, 2026
% =========================================================================

%% ------------------------------------------------------------------------
% 0. Setup
% -------------------------------------------------------------------------

% ── Which neurons to include ─────────────────────────────────────────────
% Options: all neurons, positive ramp only, negative ramp only, ramping any
% Change this mask to suit your question
neuron_mask = neuron_idx.fp_neg_ramping;   % positive ramp
% neuron_mask = abs(analysis_log.ramping_fp_flag) == 1;  % any ramp
% neuron_mask = true(height(analysis_log), 1);            % all neurons

% ── Array and timing ─────────────────────────────────────────────────────
time_vec      = -2000:2000;       % full SDF time vector
sdf_all_conds = z_fixation_array; % nNeurons x nTime x 4
% Condition indices: 2=short, 3=medium, 4=long (1=all, skip)
cond_idx      = [2, 3, 4];
cond_labels   = {'Short (725 ms)', 'Medium (1000 ms)', 'Long (1300 ms)'};
fp_durations  = [725, 1000, 1300];   % ms
fp_colors     = [0.17 0.47 0.84;     % blue  - short
                 0.11 0.69 0.48;     % green - medium
                 0.93 0.63 0.00];    % amber - long

% ── Common PCA window: last 725 ms before stim onset (exists for all FPs) 
pca_win_ms   = [-725, 0];
pca_mask     = time_vec >= pca_win_ms(1) & time_vec <= pca_win_ms(2);

% ── Full projection window per condition ─────────────────────────────────
full_win_ms  = [-1400, 200];   % show 200 ms post-stim for context
full_mask    = time_vec >= full_win_ms(1) & time_vec <= full_win_ms(2);
t_full       = time_vec(full_mask);

% ── N PCs to show (set after inspecting scree — see Section 1) ───────────
N_PCS        = 3;   % adjust after seeing scree plot

%% ------------------------------------------------------------------------
% 1. Condition-averaged population responses + PCA subspace definition
% -------------------------------------------------------------------------

% Extract condition-averaged responses for included neurons
% Shape after averaging: nNeurons x nTime (averaged across included conds
% in the PCA window, to define a condition-general subspace)

nNeurons   = length(neuron_mask);
nConds     = numel(cond_idx);
nTime_pca  = sum(pca_mask);
nTime_full = sum(full_mask);

fprintf('Defining PCA subspace:\n');
fprintf('  Neurons: %i\n', nNeurons);
fprintf('  PCA window: %i to %i ms (%i timepoints)\n', ...
    pca_win_ms(1), pca_win_ms(2), nTime_pca);

% Build matrix for PCA: (nNeurons * nConds) x nTime_pca
% Each row is one neuron's response in one condition during the PCA window
% This pools across conditions so PCA captures shared temporal structure
X_pca = nan(nNeurons * nConds, nTime_pca);

for c_i = 1:nConds
    ci       = cond_idx(c_i);
    rows     = (c_i-1)*nNeurons + (1:nNeurons);
    X_pca(rows, :) = sdf_all_conds(neuron_mask, pca_mask, ci);
end

% Remove rows with NaN
valid_rows = ~any(isnan(X_pca), 2);
X_pca_clean = X_pca(valid_rows, :);

fprintf('  Valid rows for PCA: %i of %i\n', sum(valid_rows), size(X_pca,1));

% Run PCA — operates on (observations x variables) = (neuron-conds x time)
[coeff, score, latent, ~, explained] = pca(X_pca_clean);

% Scree plot
figure('Position', [50 50 500 300], 'Renderer', 'painters');
n_show = min(20, numel(explained));
bar(1:n_show, explained(1:n_show), 'FaceColor', [0.3 0.5 0.8], 'EdgeColor','none');
hold on;
yyaxis right;
plot(1:n_show, cumsum(explained(1:n_show)), 'o-', ...
    'Color', [0.8 0.3 0.2], 'MarkerFaceColor', [0.8 0.3 0.2], ...
    'LineWidth', 1.5, 'MarkerSize', 5);
yline(80, '--', '80%', 'Color', [0.8 0.3 0.2], 'LineWidth', 1);
ylabel('Cumulative variance (%)');
yyaxis left;
ylabel('Variance explained (%)');
xlabel('PC'); title('Scree plot — set N\_PCS then re-run');
box off; set(gca, 'TickDir', 'out');

fprintf('\nVariance explained:\n');
for thresh = [70 80 90]
    n_thresh = find(cumsum(explained) >= thresh, 1);
    fprintf('  %i PCs for %i%%\n', n_thresh, thresh);
end
fprintf('\nUsing N_PCS = %i (%.1f%% variance)\n', N_PCS, ...
    sum(explained(1:N_PCS)));

% PC axes (nTime_pca x N_PCS) — the subspace we'll project into
V = coeff(:, 1:N_PCS);   % time x PC

%% ------------------------------------------------------------------------
% 2. Project full-length trajectories into PCA subspace
% -------------------------------------------------------------------------
% For each condition, take the full window population response and project
% it into the PCA space defined above.
% This gives condition trajectories of different lengths in the same space.

% traj{c_i} : nTime_full x N_PCS
traj        = cell(nConds, 1);
traj_norm   = cell(nConds, 1);   % time-normalised version
N_NORM_PTS  = 100;               % normalised time grid

for c_i = 1:nConds
    ci = cond_idx(c_i);

    % Population mean in full window: nNeurons x nTime_full
    pop_response = squeeze(mean(sdf_all_conds(neuron_mask, full_mask, ci), 1, 'omitnan'));
    % pop_response is nTime_full x 1 if nNeurons=1, fix:
    if isvector(pop_response) && nNeurons > 1
        pop_response = sdf_all_conds(neuron_mask, full_mask, ci);
        pop_response = squeeze(nanmean(pop_response, 1))';  % 1 x nTime_full
    end

    % For projection: need (nTime_full x nNeurons) * (nNeurons x N_PCS)
    % But our PCA was on time x neuron-cond, so V is nTime_pca x N_PCS
    % We need to project neuron responses onto PC time axes
    % Population trajectory = mean across neurons of (neuron response · V)

    % neuron_responses: nNeurons x nTime_full
    neuron_resp = squeeze(sdf_all_conds(neuron_mask, full_mask, ci));
    if size(neuron_resp, 1) ~= nNeurons
        neuron_resp = neuron_resp';
    end

    % Project each neuron into PC space using only the PCA window timepoints
    % then average across neurons to get population trajectory
    neuron_resp_pca = squeeze(sdf_all_conds(neuron_mask, pca_mask, ci));
    if size(neuron_resp_pca, 1) ~= nNeurons
        neuron_resp_pca = neuron_resp_pca';
    end

    % Project: for each neuron, score = response_pca_window · V
    % (nNeurons x nTime_pca) * (nTime_pca x N_PCS) = nNeurons x N_PCS
    neuron_scores = neuron_resp_pca * V;   % nNeurons x N_PCS per timepoint? No...

    % Correct approach: treat time as the variable, neurons as observations
    % Population trajectory: average neuron response at each timepoint,
    % giving 1 x nTime_full population PSTH, then project that across time
    % using V: pop_psth (nTime_full x 1) cannot project with V (nTime_pca x N_PCS)
    % Instead: project full-time population vector using V defined on PCA win
    % Standard approach: pop response is nNeurons x nTime; project as
    % traj = pop_resp_full_time' * V_neurons where V_neurons = nNeurons x N_PCS

    % Redefine: PCA on neurons x time (standard neuroscience convention)
    % See correction in Section 1b below — we redo PCA with this convention
    traj{c_i} = [];   % placeholder — filled after Section 1b
end

%% ------------------------------------------------------------------------
% 1b. Redo PCA with neuroscience convention: neurons x time
%     PCA on condition-averaged population PSTH matrix
% -------------------------------------------------------------------------
% Standard approach (Cunningham & Yu 2014, Shenoy et al.):
%   Data matrix X: nNeurons x (nTime_pca * nConds)
%   PCA gives V: nNeurons x N_PCS (neural modes)
%   Trajectory: X_full * ... projects nTime x nNeurons onto nNeurons-space
%   Actually: traj = (nTime x nNeurons) * (nNeurons x N_PCS) = nTime x N_PCS

% Build X_pca: nNeurons x (nTime_pca * nConds)
X_pca2 = nan(nNeurons, nTime_pca * nConds);
for c_i = 1:nConds
    ci   = cond_idx(c_i);
    cols = (c_i-1)*nTime_pca + (1:nTime_pca);
    data = squeeze(sdf_all_conds(neuron_mask, pca_mask, ci));
    if size(data, 1) ~= nNeurons, data = data'; end
    X_pca2(:, cols) = data;
end

% Remove neurons with NaN
valid_neurons = ~any(isnan(X_pca2), 2);
X_pca2_clean  = X_pca2(valid_neurons, :);
nNeurons_clean = sum(valid_neurons);

fprintf('\nPCA (neurons x time convention):\n');
fprintf('  Valid neurons: %i of %i\n', nNeurons_clean, nNeurons);

% PCA: observations = neurons, variables = time x condition
[coeff2, ~, latent2, ~, explained2] = pca(X_pca2_clean');
% coeff2: (nTime_pca*nConds) x N_PCS — not directly useful
% We want V_neurons: nNeurons x N_PCS (the neural axes)
% Use SVD directly for clarity:
[U, S, V_neurons_full] = svd(X_pca2_clean - mean(X_pca2_clean, 2), 'econ');
% U: nNeurons x rank, S: diagonal, V_neurons_full: (nTime*nConds) x rank
% Neural modes (columns of U): nNeurons x rank
V_neurons = U(:, 1:N_PCS);   % nNeurons x N_PCS — the subspace axes

fprintf('  Singular values accounted for: %.1f%% (first %i modes)\n', ...
    100 * sum(diag(S(1:N_PCS,1:N_PCS)).^2) / sum(diag(S).^2), N_PCS);

%% ------------------------------------------------------------------------
% 2b. Project full-length trajectories into neural subspace
% -------------------------------------------------------------------------

for c_i = 1:nConds
    ci = cond_idx(c_i);

    % Full population PSTH: nNeurons x nTime_full
    pop_full = squeeze(sdf_all_conds(neuron_mask, full_mask, ci));
    if size(pop_full, 1) ~= nNeurons, pop_full = pop_full'; end

    % Keep only valid neurons
    pop_full_clean = pop_full(valid_neurons, :);

    % Handle NaN by interpolation within each neuron
    for n_i = 1:nNeurons_clean
        row = pop_full_clean(n_i,:);
        if any(isnan(row))
            xi = find(~isnan(row));
            if numel(xi) > 1
                pop_full_clean(n_i,:) = interp1(xi, row(xi), ...
                    1:nTime_full, 'linear', 'extrap');
            else
                pop_full_clean(n_i,:) = 0;
            end
        end
    end

    % Subtract mean (same as in SVD)
    pop_full_demean = pop_full_clean - mean(X_pca2_clean, 2);

    % Project: (nTime_full x nNeurons) * (nNeurons x N_PCS)
    traj{c_i} = (pop_full_demean' * V_neurons);   % nTime_full x N_PCS

    % Time-normalised version: interpolate onto common grid
    % Define foreperiod window: -fp_duration to 0 ms
    dur       = fp_durations(c_i);
    fp_t_mask = t_full >= -dur & t_full <= 0;
    traj_fp   = traj{c_i}(fp_t_mask, :);
    t_norm    = linspace(0, 1, N_NORM_PTS);
    traj_norm{c_i} = interp1(linspace(0,1,sum(fp_t_mask)), traj_fp, t_norm);
end

%% ------------------------------------------------------------------------
% 3. Figure: trajectories aligned to stimulus onset (tests convergence)
% -------------------------------------------------------------------------

fig3 = figure('Position', [50 50 900 400], 'Renderer', 'painters');

% Panel A: PC1 vs PC2
ax1 = subplot(1, 3, 1); hold on;
title('PC1 vs PC2 — stim aligned', 'FontSize', 9);
xlabel(sprintf('PC1 (%.0f%%)', 100*diag(S(1,1))^2/sum(diag(S).^2)), 'FontSize', 8);
ylabel(sprintf('PC2'), 'FontSize', 8);
box off; set(gca, 'TickDir', 'out', 'FontSize', 8);

for c_i = 1:nConds
    plot(traj{c_i}(:,1), traj{c_i}(:,2), '-', ...
        'Color', fp_colors(c_i,:), 'LineWidth', 2, ...
        'DisplayName', cond_labels{c_i});

    % Mark fixation onset (start of foreperiod)
    dur     = fp_durations(c_i);
    fp_start = find(t_full >= -dur, 1);
    scatter(traj{c_i}(fp_start,1), traj{c_i}(fp_start,2), 40, ...
        fp_colors(c_i,:), 'filled', 'MarkerFaceAlpha', 0.5, ...
        'HandleVisibility', 'off');

    % Mark stim onset (t=0)
    t0_idx = find(t_full >= 0, 1);
    scatter(traj{c_i}(t0_idx,1), traj{c_i}(t0_idx,2), 80, ...
        fp_colors(c_i,:), 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
end
legend('Location', 'best', 'Box', 'off', 'FontSize', 7);

% Panel B: PC1 vs PC3
ax2 = subplot(1, 3, 2); hold on;
title('PC1 vs PC3 — stim aligned', 'FontSize', 9);
xlabel('PC1', 'FontSize', 8); ylabel('PC3', 'FontSize', 8);
box off; set(gca, 'TickDir', 'out', 'FontSize', 8);

for c_i = 1:nConds
    plot(traj{c_i}(:,1), traj{c_i}(:,3), '-', ...
        'Color', fp_colors(c_i,:), 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    t0_idx = find(t_full >= 0, 1);
    scatter(traj{c_i}(t0_idx,1), traj{c_i}(t0_idx,3), 80, ...
        fp_colors(c_i,:), 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
end

% Panel C: PC scores over time (each PC separately)
ax3 = subplot(1, 3, 3); hold on;
title('PC1 over time — stim aligned', 'FontSize', 9);
xlabel('Time from stim onset (ms)', 'FontSize', 8);
ylabel('PC1 projection', 'FontSize', 8);
xline(0, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
box off; set(gca, 'TickDir', 'out', 'FontSize', 8);

for c_i = 1:nConds
    plot(t_full, traj{c_i}(:,1), '-', ...
        'Color', fp_colors(c_i,:), 'LineWidth', 2, ...
        'DisplayName', cond_labels{c_i});
    % Mark foreperiod start
    dur     = fp_durations(c_i);
    fp_start_t = -dur;
    xline(fp_start_t, '--', 'Color', fp_colors(c_i,:), ...
        'LineWidth', 0.8, 'Alpha', 0.5, 'HandleVisibility', 'off');
end

sgtitle('Population trajectories — stim onset aligned', 'FontSize', 11);
exportgraphics(fig3, fullfile(dirs.figures, 'pca_traj_stim_aligned.pdf'), ...
    'ContentType', 'vector');
fprintf('Saved: pca_traj_stim_aligned.pdf\n');

%% ------------------------------------------------------------------------
% 4. Figure: time-normalised trajectories (tests temporal scaling)
% -------------------------------------------------------------------------

fig4 = figure('Position', [50 50 700 350], 'Renderer', 'painters');

subplot(1, 2, 1); hold on;
title('PC1 vs PC2 — time normalised', 'FontSize', 9);
xlabel('PC1', 'FontSize', 8); ylabel('PC2', 'FontSize', 8);
box off; set(gca, 'TickDir', 'out', 'FontSize', 8);

for c_i = 1:nConds
    plot(traj_norm{c_i}(:,1), traj_norm{c_i}(:,2), '-', ...
        'Color', fp_colors(c_i,:), 'LineWidth', 2, ...
        'DisplayName', cond_labels{c_i});
    scatter(traj_norm{c_i}(1,1),   traj_norm{c_i}(1,2),   40, ...
        fp_colors(c_i,:), 'filled', 'MarkerFaceAlpha', 0.5, ...
        'HandleVisibility', 'off');
    scatter(traj_norm{c_i}(end,1), traj_norm{c_i}(end,2), 80, ...
        fp_colors(c_i,:), 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
end
legend('Location', 'best', 'Box', 'off', 'FontSize', 7);
text(0.5, 0.02, {'If trajectories overlap:', 'temporal scaling supported'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 7, 'Color', [0.5 0.5 0.5]);

subplot(1, 2, 2); hold on;
title('PC1 over normalised time', 'FontSize', 9);
xlabel('Normalised time (0=fix, 1=stim)', 'FontSize', 8);
ylabel('PC1 projection', 'FontSize', 8);
box off; set(gca, 'TickDir', 'out', 'FontSize', 8);
t_norm_axis = linspace(0, 1, N_NORM_PTS);

for c_i = 1:nConds
    plot(t_norm_axis, traj_norm{c_i}(:,1), '-', ...
        'Color', fp_colors(c_i,:), 'LineWidth', 2, ...
        'DisplayName', cond_labels{c_i});
end
legend('Location', 'best', 'Box', 'off', 'FontSize', 7);

sgtitle('Time-normalised trajectories — temporal scaling test', 'FontSize', 11);
exportgraphics(fig4, fullfile(dirs.figures, 'pca_traj_normalised.pdf'), ...
    'ContentType', 'vector');
fprintf('Saved: pca_traj_normalised.pdf\n');

%% ------------------------------------------------------------------------
% 5. Formal test: principal angles between condition subspaces
% -------------------------------------------------------------------------
% Computes the principal angles between the low-D subspaces spanned by
% each pair of conditions in the PCA window.
% Small angles = shared subspace (consistent with common coding scheme).
% Large angles = distinct subspaces (condition-specific geometry).

fprintf('\n--- Principal angles between condition subspaces ---\n');

% For each condition, get the top N_PCS left singular vectors
% of the (nNeurons x nTime_pca) response matrix
cond_bases = cell(nConds, 1);
for c_i = 1:nConds
    ci   = cond_idx(c_i);
    data = squeeze(sdf_all_conds(neuron_mask, pca_mask, ci));
    if size(data,1) ~= nNeurons, data = data'; end
    data = data(valid_neurons, :);
    data = data - mean(data, 2);
    [U_c, ~, ~] = svd(data, 'econ');
    cond_bases{c_i} = U_c(:, 1:N_PCS);
end

pairs = nchoosek(1:nConds, 2);
for p_i = 1:size(pairs, 1)
    a = pairs(p_i, 1);
    b = pairs(p_i, 2);

    % Principal angles via SVD of cross-Gram matrix
    M      = cond_bases{a}' * cond_bases{b};
    [~,S_pa,~] = svd(M);
    angles_deg = acosd(min(diag(S_pa), 1));

    fprintf('  %s vs %s:\n', cond_labels{a}, cond_labels{b});
    fprintf('    Principal angles: %s degrees\n', ...
        strjoin(arrayfun(@(x) sprintf('%.1f', x), angles_deg', ...
        'UniformOutput', false), ', '));
    fprintf('    Mean angle: %.1f degrees\n', mean(angles_deg));
    fprintf('    (0° = identical subspace, 90° = orthogonal)\n');
end