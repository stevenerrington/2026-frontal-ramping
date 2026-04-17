
% Extract activity aligned on stimulus onset, across multiple foreperiods
fixation_array = [];

for neuron_i = 1:size(spike_log,1)

    fprintf('Analysing neuron %i of %i \n', neuron_i, size(spike_log,1));
    neuron_label = spike_log.neuron_label{neuron_i};
    session_idx = find(strcmp(session_log.session, spike_log.session{neuron_i}));

    beh_in = load(fullfile(dirs.raw_data, [spike_log.session{neuron_i} ,'_spk.mat']));
    sdf_in = load(fullfile('/Volumes/Mnemosyne/Data/2026_macaque_value/proc/', [neuron_label '.mat']));

    clear ttx foreperiod
    foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;
    ttx.short_fp = find(foreperiod > 0.7 & foreperiod < 0.75);
    ttx.mid_fp = find(foreperiod > 0.9 & foreperiod < 1.10);
    ttx.long_fp = find(foreperiod > 1.2);
    ttx.all_fp = find(foreperiod > 0.5);


    fixation_array(neuron_i,:,1) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.all_fp ,:)),100)';
    fixation_array(neuron_i,:,2) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.short_fp ,:)),100)';
    fixation_array(neuron_i,:,3) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.mid_fp ,:)),100)';
    fixation_array(neuron_i,:,4) = smooth(nanmean(sdf_in.sdf.stim_on(ttx.long_fp ,:)),100)';

end

save(fullfile(dirs.root,'_data','sdf','fixation_array.mat'), 'fixation_array', '-v7.3')


%% 

linearity_out_r2 = []; linearity_out_pval = [];
parfor neuron_i = 1:size(spike_log,1)
    disp(['Assessing linearity in neuron ' int2str(neuron_i)])
    sdf_zero = 2000;
    ramp_win = [-750:0];
    ramp_sdf = fixation_array(neuron_i,sdf_zero+ramp_win,1);
    
    null_R2 = [];
    for i = 1:nIter
        sdf_null = circshift(ramp_sdf, randi(length(ramp_sdf)));
        mdl = fitlm(t, sdf_null);
        null_R2(i) = mdl.Rsquared.Ordinary;
    end

    mdl_real = fitlm(ramp_win, ramp_sdf);
    real_R2 = mdl_real.Rsquared.Ordinary;
    pval = mean(null_R2 >= real_R2);

    linearity_out_r2(neuron_i,1) = real_R2;
    linearity_out_pval(neuron_i, 2) = pval;
end