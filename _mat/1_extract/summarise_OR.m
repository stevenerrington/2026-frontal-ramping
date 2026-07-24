%% ========================================================================
% extract_or_results.m
% Extracts all statistics needed for the odds ratio results section.
% Run this after ramp_prevalence_OR.m — requires the same workspace.
% Prints a structured summary ready to paste into a results write-up.
% Steven Errington, 2026
% =========================================================================

fprintf('\n========================================================\n');
fprintf(' ODDS RATIO ANALYSIS — RESULTS SUMMARY\n');
fprintf('========================================================\n\n');

epochs = struct( ...
    'label', {'Foreperiod',      'Reward'}, ...
    'flag',  {'ramping_fp_flag', 'ramping_rwd_flag'});

nEpochs = numel(epochs);

for ep_i = 1:nEpochs

    flag_col  = epochs(ep_i).flag;
    ep_label  = epochs(ep_i).label;
    flags     = analysis_log.(flag_col);

    n_total   = height(analysis_log);
    n_pos     = sum(flags ==  1);
    n_neg     = sum(flags == -1);
    n_ramp    = n_pos + n_neg;
    n_nonramp = sum(flags ==  0);

    fprintf('--------------------------------------------------------\n');
    fprintf('EPOCH: %s\n', ep_label);
    fprintf('--------------------------------------------------------\n');
    fprintf('Total neurons     : %i\n',         n_total);
    fprintf('Ramping (any)     : %i (%.1f%%)\n', n_ramp,    100*n_ramp/n_total);
    fprintf('  Positive ramp   : %i (%.1f%%)\n', n_pos,     100*n_pos/n_total);
    fprintf('  Negative ramp   : %i (%.1f%%)\n', n_neg,     100*n_neg/n_total);
    fprintf('Non-ramping       : %i (%.1f%%)\n', n_nonramp, 100*n_nonramp/n_total);

    % ── Omnibus chi-square across areas ──────────────────────────────────
    cont_area = zeros(nPlot, 2);
    for area_i = 1:nPlot
        area_mask = strcmp(analysis_log.area, area_order_plot{area_i});
        cont_area(area_i, 1) = sum(abs(flags(area_mask)) == 1);   % ramp
        cont_area(area_i, 2) = sum(flags(area_mask) == 0);         % nonramp
    end
    row_s    = sum(cont_area, 2);
    col_s    = sum(cont_area, 1);
    N        = sum(cont_area(:));
    expected = (row_s * col_s) / N;
    chi2_omni = sum((cont_area(:) - expected(:)).^2 ./ expected(:));
    df_omni   = (size(cont_area,1)-1) * (size(cont_area,2)-1);
    p_omni    = 1 - chi2cdf(chi2_omni, df_omni);

    fprintf('\nOmnibus chi-square (ramp prevalence across areas):\n');
    fprintf('  X²(%i) = %.2f, p = %.2e\n', df_omni, chi2_omni, p_omni);

    % ── Per-area OR results ───────────────────────────────────────────────
    % Pull from the pre-computed matrices (from ramp_prevalence_OR.m)
    ep_or     = or_val(:, ep_i);
    ep_ci_lo  = ci_lo(:, ep_i);
    ep_ci_hi  = ci_hi(:, ep_i);
    ep_p_fdr  = p_fdr(:, ep_i);
    ep_p_raw  = chi2_p(:, ep_i);

    sig_mask  = ep_p_fdr < 0.05;
    n_sig     = sum(sig_mask & ~isnan(ep_p_fdr));
    n_tested  = sum(~isnan(ep_p_fdr));

    fprintf('\nPer-area OR (one-vs-rest, FDR corrected):\n');
    fprintf('  Significant areas: %i of %i tested\n', n_sig, n_tested);

    % Enriched (OR > 1, significant)
    enriched = find(sig_mask & ep_or > 1);
    if ~isempty(enriched)
        fprintf('  Enriched (OR > 1):\n');
        for i = 1:numel(enriched)
            a = enriched(i);
            fprintf('    %-8s  OR = %.2f [%.2f–%.2f]  p_raw = %.4f  p_fdr = %.4f\n', ...
                area_order_plot{a}, ep_or(a), ep_ci_lo(a), ep_ci_hi(a), ...
                ep_p_raw(a), ep_p_fdr(a));
        end
    end

    % Depleted (OR < 1, significant)
    depleted = find(sig_mask & ep_or < 1);
    if ~isempty(depleted)
        fprintf('  Depleted (OR < 1):\n');
        for i = 1:numel(depleted)
            a = depleted(i);
            fprintf('    %-8s  OR = %.2f [%.2f–%.2f]  p_raw = %.4f  p_fdr = %.4f\n', ...
                area_order_plot{a}, ep_or(a), ep_ci_lo(a), ep_ci_hi(a), ...
                ep_p_raw(a), ep_p_fdr(a));
        end
    end

    % Non-significant areas (print more briefly)
    nonsig = find(~sig_mask & ~isnan(ep_p_fdr));
    if ~isempty(nonsig)
        fprintf('  Non-significant areas: ');
        fprintf('%s ', area_order_plot{nonsig});
        fprintf('\n');
    end

    fprintf('\n');
end

% ── Overlap between epochs ────────────────────────────────────────────────
fprintf('--------------------------------------------------------\n');
fprintf('EPOCH OVERLAP (collapsed across areas)\n');
fprintf('--------------------------------------------------------\n');

ramp_fp  = abs(analysis_log.ramping_fp_flag)  == 1;
ramp_rwd = abs(analysis_log.ramping_rwd_flag) == 1;

n_total    = height(analysis_log);
n_fp_only  = sum( ramp_fp & ~ramp_rwd);
n_rwd_only = sum(~ramp_fp &  ramp_rwd);
n_both     = sum( ramp_fp &  ramp_rwd);
n_neither  = sum(~ramp_fp & ~ramp_rwd);

fprintf('Foreperiod only    : %i (%.1f%%)\n', n_fp_only,  100*n_fp_only/n_total);
fprintf('Reward only        : %i (%.1f%%)\n', n_rwd_only, 100*n_rwd_only/n_total);
fprintf('Both epochs        : %i (%.1f%%)\n', n_both,     100*n_both/n_total);
fprintf('Neither            : %i (%.1f%%)\n', n_neither,  100*n_neither/n_total);

cont_overlap = [n_both, n_fp_only; n_rwd_only, n_neither];
[~, p_fisher] = fishertest(cont_overlap);
fprintf('Fisher exact test  : p = %.4f\n', p_fisher);

% Expected overlap under independence
p_fp  = (n_fp_only + n_both) / n_total;
p_rwd = (n_rwd_only + n_both) / n_total;
n_expected_both = n_total * p_fp * p_rwd;
fprintf('Expected overlap (independence): %.0f neurons\n', n_expected_both);
fprintf('Observed overlap               : %i neurons\n', n_both);
fprintf('Enrichment ratio               : %.2fx\n', n_both / n_expected_both);

fprintf('\n========================================================\n');
fprintf(' END OF SUMMARY\n');
fprintf('========================================================\n\n');