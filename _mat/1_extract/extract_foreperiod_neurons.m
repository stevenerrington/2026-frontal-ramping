%% ========================================================================
% Extract SDF activity aligned to stimulus onset, split by foreperiod
% =========================================================================

savePath = fullfile(dirs.root, '_data', 'sdf', 'fixation_array.mat');
nNeurons = size(spike_log, 1);

if exist(savePath, 'file')
    fprintf('Loading existing fixation_array from disk...\n')
    load(savePath, 'fixation_array')
else
    fprintf('Computing fixation_array for %i neurons...\n', nNeurons)
    fixation_array = [];

    for neuron_i = 1:nNeurons
        fprintf('  Neuron %i of %i\n', neuron_i, nNeurons)

        neuron_label = spike_log.neuron_label{neuron_i};

        % Load behavioural and SDF data
        beh_in = load(fullfile(dirs.raw_data, [spike_log.session{neuron_i}, '_spk.mat']));
        sdf_in = load(fullfile(dirs.proc_data, [neuron_label, '.mat']));

        % Compute foreperiod duration (fixation cross fix to stim on)
        foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;

        % Trial indices by foreperiod duration
        ttx.all_fp   = foreperiod > 0.5;
        ttx.short_fp = foreperiod > 0.70 & foreperiod < 0.75;
        ttx.mid_fp   = foreperiod > 0.90 & foreperiod < 1.10;
        ttx.long_fp  = foreperiod > 1.20;

        % Smooth trial-averaged SDF for each foreperiod condition
        % Dim 3: 1 = all, 2 = short, 3 = medium, 4 = long
        fixation_array(neuron_i, :, 1) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.all_fp,   :), 1), 100)';
        fixation_array(neuron_i, :, 2) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.short_fp, :), 1), 100)';
        fixation_array(neuron_i, :, 3) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.mid_fp,   :), 1), 100)';
        fixation_array(neuron_i, :, 4) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.long_fp,  :), 1), 100)';
    end

    save(savePath, 'fixation_array', '-v7.3')
    fprintf('fixation_array saved to %s\n', savePath)
end

%% ========================================================================
% Extract SDF activity aligned to fixation onset, split by foreperiod
% =========================================================================

savePath = fullfile(dirs.root, '_data', 'sdf', 'fixon_array.mat');
nNeurons = size(spike_log, 1);

if exist(savePath, 'file')
    fprintf('Loading existing fixon_array from disk...\n')
    load(savePath, 'fixon_array')
else
    fprintf('Computing fixon_array for %i neurons...\n', nNeurons)
    fixon_array = [];

    for neuron_i = 1:nNeurons
        fprintf('  Neuron %i of %i\n', neuron_i, nNeurons)

        neuron_label = spike_log.neuron_label{neuron_i};

        % Load behavioural and SDF data
        beh_in = load(fullfile(dirs.raw_data, [spike_log.session{neuron_i}, '_spk.mat']));
        sdf_in = load(fullfile(dirs.proc_data, [neuron_label, '.mat']));

        % Compute foreperiod duration (fixation cross fix to stim on)
        foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;

        % Trial indices by foreperiod duration
        ttx.all_fp   = foreperiod > 0.5;
        ttx.short_fp = foreperiod > 0.70 & foreperiod < 0.75;
        ttx.mid_fp   = foreperiod > 0.90 & foreperiod < 1.10;
        ttx.long_fp  = foreperiod > 1.20;

        % Smooth trial-averaged SDF for each foreperiod condition
        % Dim 3: 1 = all, 2 = short, 3 = medium, 4 = long
        fixon_array(neuron_i, :, 1) = smooth(nanmean(sdf_in.sdf.fixation(ttx.all_fp,   :), 1), 100)';
        fixon_array(neuron_i, :, 2) = smooth(nanmean(sdf_in.sdf.fixation(ttx.short_fp, :), 1), 100)';
        fixon_array(neuron_i, :, 3) = smooth(nanmean(sdf_in.sdf.fixation(ttx.mid_fp,   :), 1), 100)';
        fixon_array(neuron_i, :, 4) = smooth(nanmean(sdf_in.sdf.fixation(ttx.long_fp,  :), 1), 100)';
    end

    save(savePath, 'fixon_array', '-v7.3')
    fprintf('fixon_array saved to %s\n', savePath)
end

%% ========================================================================
% Z-score SDF for each neuron and condition
% =========================================================================

nConds        = size(fixation_array, 3);
z_fixation_array = zeros(size(fixation_array));

for neuron_i = 1:length(fixation_array)
    z_fixation_array(neuron_i, :, 1) = zscore(fixation_array(neuron_i, :, 1));

    split_cond_array = [];
    split_cond_array = [fixation_array(neuron_i, :, 2); fixation_array(neuron_i, :, 3); fixation_array(neuron_i, :, 4)];

    z_fixation_array(neuron_i, :, 2) =  (fixation_array(neuron_i, :, 2) - mean(split_cond_array(:)))./std(split_cond_array(:));
    z_fixation_array(neuron_i, :, 3) =  (fixation_array(neuron_i, :, 3) - mean(split_cond_array(:)))./std(split_cond_array(:));
    z_fixation_array(neuron_i, :, 4) =  (fixation_array(neuron_i, :, 4) - mean(split_cond_array(:)))./std(split_cond_array(:));

end

%% ========================================================================
% Assess foreperiod ramping via linear fit + circular shift permutation
% =========================================================================

fp_linearity_file = fullfile(dirs.root, '_data', 'fp_linearity_analysis.csv');

if exist(fp_linearity_file, 'file')
    fprintf('Loading datatable...\n')
    fp_linearity_table = readtable(fp_linearity_file);
else
    sdf_zero  = 2000;                  % index of time zero (stim_on) in SDF
    ramp_ms   = -750:0;               % ramp window in ms relative to stim_on
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

    linearity_analysis_table = array2table(linearity_out, 'VariableNames', {'flag','pval','slope','r2'});

    linearity_analysis_table.area = spike_log.area;
    linearity_analysis_table.neuron_label = spike_log.neuron_label;

    writetable(linearity_analysis_table, fp_linearity_file, 'WriteRowNames', true);
end

% ------------------------------------------------------------------------
% Classify ramping neurons and append flag to analysis_log
% -------------------------------------------------------------------------

neuron_idx.fp_pos_ramping = find(...
    fp_linearity_table.flag  == 1  & ...
    fp_linearity_table.slope >  0  & ...
    fp_linearity_table.r2    >  0.8);

neuron_idx.fp_neg_ramping = find(...
    fp_linearity_table.flag  == 1  & ...
    fp_linearity_table.slope <  0  & ...
    fp_linearity_table.r2    >  0.8);

fprintf('Positive ramp neurons: %i\n', length(neuron_idx.fp_pos_ramping))
fprintf('Negative ramp neurons: %i\n', length(neuron_idx.fp_neg_ramping))

% Append ramping flag to log (1 = pos, -1 = neg, 0 = non-ramping)
analysis_log.ramping_fp_flag = zeros(size(analysis_log, 1), 1);
analysis_log.ramping_fp_flag(neuron_idx.fp_pos_ramping) =  1;
analysis_log.ramping_fp_flag(neuron_idx.fp_neg_ramping) = -1;