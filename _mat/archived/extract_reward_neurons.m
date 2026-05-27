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

        valid_trials = find(beh_in.trialtype.rew == 1);

        reward_array(neuron_i, :) = smooth(nanmean(sdf_in.sdf.stim_on(valid_trials,   :), 1), 100)';
    end

    save(savePath, 'reward_array', '-v7.3')
    fprintf('reward_array saved to %s\n', savePath)
end

%% ========================================================================
% Z-score SDF for each neuron and condition
% =========================================================================

nConds        = size(reward_array, 3);
z_fixation_array = zeros(size(reward_array));

for neuron_i = 1:nNeurons
        z_fixation_array(neuron_i, :) = zscore(reward_array(neuron_i, :));
end

%% ========================================================================
% Assess foreperiod ramping via linear fit + circular shift permutation
% =========================================================================

sdf_zero  = 2000;                  % index of time zero (stim_on) in SDF
ramp_ms   = -400:0;               % ramp window in ms relative to stim_on
ramp_idx  = sdf_zero + ramp_ms;   % corresponding SDF indices
t         = ramp_ms;              % time vector for fitting
nPerm     = 1000;
nPts      = length(t);

% Pre-allocate: cols = [isRamping, pVal, slope, Rsq]
linearity_out = nan(nNeurons, 4);

for neuron_i = 1:nNeurons
    fprintf('  Assessing linearity: neuron %i of %i\n', neuron_i, nNeurons)

    fr = z_fixation_array(neuron_i, ramp_idx, 1);  % all-FP z-scored SDF

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

reward_linearity_analysis_table = array2table(linearity_out, 'VariableNames', {'flag','pval','slope','r2'});

reward_linearity_analysis_table.area = spike_log.area;
reward_linearity_analysis_table.neuron_label = spike_log.neuron_label;

writetable(reward_linearity_analysis_table, fullfile(dirs.root,'_data','reward_linearity_analysis_table.csv'), 'WriteRowNames', true);