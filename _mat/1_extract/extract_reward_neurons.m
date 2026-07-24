%% ========================================================================
% Extract SDF activity aligned to stimulus onset, split by foreperiod
% =========================================================================

savePath = fullfile(dirs.root, '_data', 'sdf', 'reward_array.mat');
nNeurons = size(spike_log, 1);

if exist(savePath, 'file')
    fprintf('Loading existing reward_array from disk...\n')
    load(savePath, 'reward_array')
else
    fprintf('Computing reward_array for %i neurons...\n', nNeurons)
    reward_array = [];

    for neuron_i = 1:nNeurons
        fprintf('  Neuron %i of %i\n', neuron_i, nNeurons)

        neuron_label = spike_log.neuron_label{neuron_i};

        % Load behavioural and SDF data
        beh_in = load(fullfile(dirs.raw_data, [spike_log.session{neuron_i}, '_spk.mat']));
        sdf_in = load(fullfile(dirs.proc_data, [neuron_label, '.mat']));

        valid_trials = beh_in.trialtype.rew == 1;
        low_prob = find(valid_trials & beh_in.trialtype.chosenproba_2AFC < 35);
        mid_prob = find(valid_trials & beh_in.trialtype.chosenproba_2AFC == 50);
        high_prob = find(valid_trials & beh_in.trialtype.chosenproba_2AFC > 65);

        reward_array(neuron_i, :, 1) = smooth(nanmean(sdf_in.sdf.reward(valid_trials,   :), 1), 100)';
        reward_array(neuron_i, :, 2) = smooth(nanmean(sdf_in.sdf.reward(low_prob,   :), 1), 100)';
        reward_array(neuron_i, :, 3) = smooth(nanmean(sdf_in.sdf.reward(mid_prob,   :), 1), 100)';
        reward_array(neuron_i, :, 4) = smooth(nanmean(sdf_in.sdf.reward(high_prob,   :), 1), 100)';
    end

    save(savePath, 'reward_array', '-v7.3')
    fprintf('reward_array saved to %s\n', savePath)
end

%% ========================================================================
% Z-score SDF for each neuron and condition
% =========================================================================

nConds        = size(reward_array, 3);
z_reward_array = zeros(size(reward_array));

for neuron_i = 1:length(reward_array)
    z_reward_array(neuron_i, :, 1) = zscore(reward_array(neuron_i, :, 1));

    split_cond_array = [];
    split_cond_array = [reward_array(neuron_i, :, 2); reward_array(neuron_i, :, 3); reward_array(neuron_i, :, 4)];

    z_reward_array(neuron_i, :, 2) =  (reward_array(neuron_i, :, 2) - mean(split_cond_array(:)))./std(split_cond_array(:));
    z_reward_array(neuron_i, :, 3) =  (reward_array(neuron_i, :, 3) - mean(split_cond_array(:)))./std(split_cond_array(:));
    z_reward_array(neuron_i, :, 4) =  (reward_array(neuron_i, :, 4) - mean(split_cond_array(:)))./std(split_cond_array(:));

end

%% ========================================================================
% Assess foreperiod ramping via linear fit + circular shift permutation
% =========================================================================

reward_linearity_file = fullfile(dirs.root, '_data', 'reward_linearity_analysis.csv');

if exist(reward_linearity_file, 'file')
    fprintf('Loading datatable...\n')
    reward_linearity_table = readtable(reward_linearity_file);
else
    sdf_zero  = 2000;                  % index of time zero (stim_on) in SDF
    ramp_ms   = -1000:0;               % ramp window in ms relative to stim_on
    ramp_idx  = sdf_zero + ramp_ms;   % corresponding SDF indices
    t         = ramp_ms;              % time vector for fitting
    nPerm     = 1000;
    nPts      = length(t);

    % Pre-allocate: cols = [isRamping, pVal, slope, Rsq]
    linearity_out = nan(nNeurons, 4);

    for neuron_i = 1:nNeurons
        fprintf('  Assessing linearity: neuron %i of %i\n', neuron_i, nNeurons)

        fr = z_reward_array(neuron_i, ramp_idx, 1);  % all-FP z-scored SDF

        % Linear fit
        p    = polyfit(t, fr, 1);
        slope = p(1);
        yhat  = polyval(p, t);
        ss_res = sum((fr - yhat).^2);
        ss_tot = sum((fr - mean(fr)).^2);
        Rsq    = 1 - ss_res / ss_tot;

        % Circular shift permutation null distribution
        RsqNull = zeros(nPerm, 1);
        for perm_i = 1:nPerm
            shiftAmt      = randi(nPts - 1);
            frShift       = circshift(fr, shiftAmt);
            pNull         = polyfit(t, frShift, 1);
            yhatNull      = polyval(pNull, t);
            ss_resNull    = sum((frShift - yhatNull).^2);
            ss_totNull    = sum((frShift - mean(frShift)).^2);
            RsqNull(perm_i) = 1 - ss_resNull / ss_totNull;
        end

        pVal      = mean(RsqNull >= Rsq);
        isRamping = pVal < 0.001;

        linearity_out(neuron_i, :) = [double(isRamping), pVal, slope, Rsq];
    end

    reward_linearity_table = array2table(linearity_out, 'VariableNames', {'flag','pval','slope','r2'});

    reward_linearity_table.area = spike_log.area;
    reward_linearity_table.neuron_label = spike_log.neuron_label;

    writetable(reward_linearity_table, fullfile(dirs.root,'_data','reward_linearity_analysis.csv'), 'WriteRowNames', true);
end

% ------------------------------------------------------------------------
% Classify ramping neurons and append flag to analysis_log
% -------------------------------------------------------------------------

neuron_idx.rwd_pos_ramping = find(...
    reward_linearity_table.flag  == 1  & ...
    reward_linearity_table.slope >  0  & ...
    reward_linearity_table.r2    >  0.8);

neuron_idx.rwd_neg_ramping = find(...
    reward_linearity_table.flag  == 1  & ...
    reward_linearity_table.slope <  0  & ...
    reward_linearity_table.r2    >  0.8);

fprintf('Positive ramp neurons: %i\n', length(neuron_idx.rwd_pos_ramping))
fprintf('Negative ramp neurons: %i\n', length(neuron_idx.rwd_neg_ramping))

% Append ramping flag to log (1 = pos, -1 = neg, 0 = non-ramping)
analysis_log.ramping_rwd_flag = zeros(size(analysis_log, 1), 1);
analysis_log.ramping_rwd_flag(neuron_idx.rwd_pos_ramping) =  1;
analysis_log.ramping_rwd_flag(neuron_idx.rwd_neg_ramping) = -1;

