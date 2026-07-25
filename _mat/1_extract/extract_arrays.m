%% ========================================================================
% extract_arrays.m
%
% Preprocessing step: extract the three condition-split, trial-averaged
% SDF arrays that serve as the foundation for all downstream analysis.
% This script performs EXTRACTION ONLY — no z-scoring, ramp detection,
% classification, or plotting. All of that now happens in Python
% (see frontal_ramping_analysis.ipynb).
%
% For each neuron, behavioural and SDF data are loaded once and used to
% populate all three arrays in a single pass:
%
%   fixation_array.mat  — sdf.stim_on,  split by foreperiod duration
%   fixon_array.mat     — sdf.fixation, split by foreperiod duration
%   reward_array.mat    — sdf.reward,   split by chosen reward probability
%
% Array shape: nNeurons x nTimepoints x 4
%   fixation_array / fixon_array, dim 3: 1 = all, 2 = short, 3 = medium, 4 = long
%   reward_array,                 dim 3: 1 = all, 2 = low,   3 = mid,    4 = high
%
% All three are saved with '-v7.3' for consistency and to support
% files > 2GB.
%
% Steven Errington, 2026
% =========================================================================

nNeurons = size(spike_log, 1);

fixation_path = fullfile(dirs.root, '_data', 'sdf', 'fixation_array.mat');
fixon_path    = fullfile(dirs.root, '_data', 'sdf', 'fixon_array.mat');
reward_path   = fullfile(dirs.root, '_data', 'sdf', 'reward_array.mat');

already_done = exist(fixation_path, 'file') && ...
               exist(fixon_path,    'file') && ...
               exist(reward_path,   'file');

if already_done
    fprintf('fixation_array, fixon_array, and reward_array already exist on disk.\n')
    fprintf('Delete them from _data/sdf/ and rerun this script to recompute.\n')
    return
end

fprintf('Extracting fixation_array, fixon_array, and reward_array for %i neurons...\n', nNeurons)

fixation_array = nan(nNeurons, 4001, 4);
fixon_array    = nan(nNeurons, 4001, 4);
reward_array   = nan(nNeurons, 4001, 4);

for neuron_i = 1:nNeurons
    fprintf('  Neuron %i of %i\n', neuron_i, nNeurons)

    neuron_label = spike_log.neuron_label{neuron_i};

    % Load behavioural and SDF data once per neuron
    beh_in = load(fullfile(dirs.raw_data, [spike_log.session{neuron_i}, '_spk.mat']));
    sdf_in = load(fullfile(dirs.proc_data, [neuron_label, '.mat']));

    % ── Foreperiod trial masks (shared by fixation_array & fixon_array) ────
    foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;

    ttx.all_fp   = foreperiod > 0.5;
    ttx.short_fp = foreperiod > 0.70 & foreperiod < 0.75;
    ttx.mid_fp   = foreperiod > 0.90 & foreperiod < 1.10;
    ttx.long_fp  = foreperiod > 1.20;

    % Dim 3: 1 = all, 2 = short, 3 = medium, 4 = long
    fixation_array(neuron_i, :, 1) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.all_fp,   :), 1), 100)';
    fixation_array(neuron_i, :, 2) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.short_fp, :), 1), 100)';
    fixation_array(neuron_i, :, 3) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.mid_fp,   :), 1), 100)';
    fixation_array(neuron_i, :, 4) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.long_fp,  :), 1), 100)';

    fixon_array(neuron_i, :, 1) = smooth(nanmean(sdf_in.sdf.fixation(ttx.all_fp,   :), 1), 100)';
    fixon_array(neuron_i, :, 2) = smooth(nanmean(sdf_in.sdf.fixation(ttx.short_fp, :), 1), 100)';
    fixon_array(neuron_i, :, 3) = smooth(nanmean(sdf_in.sdf.fixation(ttx.mid_fp,   :), 1), 100)';
    fixon_array(neuron_i, :, 4) = smooth(nanmean(sdf_in.sdf.fixation(ttx.long_fp,  :), 1), 100)';

    % ── Reward trial masks ──────────────────────────────────────────────────
    valid_trials = beh_in.trialtype.rew == 1;
    low_prob     = valid_trials & beh_in.trialtype.chosenproba_2AFC < 35;
    mid_prob     = valid_trials & beh_in.trialtype.chosenproba_2AFC == 50;
    high_prob    = valid_trials & beh_in.trialtype.chosenproba_2AFC > 65;

    % Dim 3: 1 = all, 2 = low, 3 = mid, 4 = high
    reward_array(neuron_i, :, 1) = smooth(nanmean(sdf_in.sdf.reward(valid_trials, :), 1), 100)';
    reward_array(neuron_i, :, 2) = smooth(nanmean(sdf_in.sdf.reward(low_prob,     :), 1), 100)';
    reward_array(neuron_i, :, 3) = smooth(nanmean(sdf_in.sdf.reward(mid_prob,     :), 1), 100)';
    reward_array(neuron_i, :, 4) = smooth(nanmean(sdf_in.sdf.reward(high_prob,    :), 1), 100)';
end

save(fixation_path, 'fixation_array', '-v7.3')
fprintf('fixation_array saved to %s\n', fixation_path)

save(fixon_path, 'fixon_array', '-v7.3')
fprintf('fixon_array saved to %s\n', fixon_path)

save(reward_path, 'reward_array', '-v7.3')
fprintf('reward_array saved to %s\n', reward_path)

clear ttx foreperiod valid_trials low_prob mid_prob high_prob beh_in sdf_in
