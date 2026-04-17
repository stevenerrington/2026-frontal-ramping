
for neuron_i = 1:20;

neuron_label = spike_log.neuron_label{neuron_i};
session_idx = find(strcmp(session_log.session, spike_log.session{neuron_i}));

beh_in = load(fullfile(dirs.raw_data, [spike_log.session{neuron_i} ,'_spk.mat']));
sdf_in = load(fullfile('/Volumes/Mnemosyne/Data/2026_macaque_value/proc/', [neuron_label '.mat']));


foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;

ttx.short_fp = find(foreperiod > 0.7 & foreperiod < 0.75);
ttx.mid_fp = find(foreperiod > 0.9 & foreperiod < 1.10);
ttx.long_fp = find(foreperiod > 1.2);



figure
subplot(1,3,1); hold on
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.fixation(ttx.short_fp,:)), 50))
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.fixation(ttx.mid_fp,:)), 50))
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.fixation(ttx.long_fp,:)), 50))
xlim([-200 1500])
vline(0, 'k-')

subplot(1,3,2); hold on
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.stim_on(:,:)), 50), 'k', 'LineWidth', 2)
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.stim_on(ttx.short_fp,:)), 50))
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.stim_on(ttx.mid_fp,:)), 50))
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.stim_on(ttx.long_fp,:)), 50))
xlim([-1500 200])
vline([-1300 -1000 -700]); vline(0, 'k-')

subplot(1,3,3)
plot(-2000:2000, smooth(nanmean(sdf_in.sdf.reward(:,:)), 50))
xlim([-600 200]); vline(0, 'k-')

end