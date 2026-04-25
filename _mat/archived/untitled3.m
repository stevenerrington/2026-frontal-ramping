
pos_neuron_i = neuron_idx.pos_ramping(35);
neg_neuron_i = neuron_idx.neg_ramping(403);



%% Figure: split x foreperiod
clear single_unit_fig

example_neurons = [pos_neuron_i, neg_neuron_i];
xlim_vals = [-1500 250];
ylim_vals = {[0 20], [0 30]};
color_pal.foreperiod_conds = slanCM('Blues',100);
color_pal.foreperiod_conds = color_pal.foreperiod_conds([33, 66, 99],:);


for example_neuron_i = 1:2

    neuron_label = spike_log.neuron_label{example_neurons(example_neuron_i)};

    % Load behavioural and SDF data
    beh_in = load(fullfile(dirs.raw_data, [spike_log.session{example_neurons(example_neuron_i)}, '_spk.mat']));
    sdf_in = load(fullfile(dirs.proc_data, [neuron_label, '.mat']));

    % Smooth SDF
    for trial_i = 1:size(sdf_in.sdf.stim_on,1)
        sdf_in.sdf.stim_on(trial_i,:) = smooth(sdf_in.sdf.stim_on(trial_i,:), 100)';
    end

    % Compute foreperiod duration (fixation cross fix to stim on)
    foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;

    % Trial indices by foreperiod duration
    ttx.all_fp   = find(foreperiod > 0.5);
    ttx.short_fp = find(foreperiod > 0.70 & foreperiod < 0.75);
    ttx.mid_fp   = find(foreperiod > 0.90 & foreperiod < 1.10);
    ttx.long_fp  = find(foreperiod > 1.20);

    trial_labels = [repmat({'Short'}, length(ttx.short_fp), 1); repmat({'Medium'}, length(ttx.mid_fp), 1); repmat({'Long'}, length(ttx.long_fp), 1)];

    % Raster plot
    single_unit_fig(1, example_neuron_i) = gramm('x', sdf_in.raster.stim_on([ttx.short_fp; ttx.mid_fp; ttx.long_fp]), 'color', trial_labels);
    single_unit_fig(1, example_neuron_i).geom_raster('geom', {'point'});
    single_unit_fig(1, example_neuron_i).axe_property('XLim', xlim_vals);
    single_unit_fig(1, example_neuron_i).set_color_options('map', color_pal.foreperiod_conds);

    % Example neuron SDF
    single_unit_fig(2, example_neuron_i) = gramm('x', -2000:2000, 'y', sdf_in.sdf.stim_on([ttx.short_fp; ttx.mid_fp; ttx.long_fp],:), 'color', trial_labels);
    single_unit_fig(2, example_neuron_i).stat_summary();
    single_unit_fig(2, example_neuron_i).axe_property('XLim', xlim_vals, 'YLim', ylim_vals{example_neuron_i});
    single_unit_fig(2, example_neuron_i).set_color_options('map', color_pal.foreperiod_conds);

end

% Add timing lines
single_unit_fig(2,1).geom_vline('xintercept',[-1300 -1000 -700 0]);
single_unit_fig(2,2).geom_vline('xintercept',[-1300 -1000 -700 0]);

% Figure setup    ////////////////////////////////////////////////////////
single_unit_fig(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
single_unit_fig(1,2).axe_property('XTick',[],'XColor',[1 1 1]);

single_unit_fig(1,1).set_layout_options...
    ('Position',[0.07 0.7 0.4 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(1,2).set_layout_options...
    ('Position',[0.55 0.7 0.4 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,1).set_layout_options...
    ('Position',[0.07 0.1 0.4 0.5],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,2).set_layout_options...
    ('Position',[0.55 0.1 0.4 0.5],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


figure('Renderer', 'painters', 'Position', [100 100 600 500]);
single_unit_fig.draw();

saveas(gcf,fullfile(dirs.figures, 'example_foreperiod_split_sdf.pdf'))


%% Figure: all foreperiods
clear single_unit_fig

xlim_vals = [-1500 250];
ylim_vals = {[0 20], [0 30]};
color_pal.foreperiod_conds = slanCM('Blues',100);
color_pal.foreperiod_conds = color_pal.foreperiod_conds([33, 66, 99],:);


for example_neuron_i = 1:2

    neuron_label = spike_log.neuron_label{example_neurons(example_neuron_i)};

    % Load behavioural and SDF data
    beh_in = load(fullfile(dirs.raw_data, [spike_log.session{example_neurons(example_neuron_i)}, '_spk.mat']));
    sdf_in = load(fullfile(dirs.proc_data, [neuron_label, '.mat']));

    % Smooth SDF
    for trial_i = 1:size(sdf_in.sdf.stim_on,1)
        sdf_in.sdf.stim_on(trial_i,:) = smooth(sdf_in.sdf.stim_on(trial_i,:), 100)';
    end

    % Compute foreperiod duration (fixation cross fix to stim on)
    foreperiod = beh_in.t_evt.stim_on - beh_in.t_evt.fixcross_fix;

    % Trial indices by foreperiod duration
    ttx.all_fp   = find(foreperiod > 0.5);
    ttx.short_fp = find(foreperiod > 0.70 & foreperiod < 0.75);
    ttx.mid_fp   = find(foreperiod > 0.90 & foreperiod < 1.10);
    ttx.long_fp  = find(foreperiod > 1.20);

    trial_labels = [repmat({'All'}, length(ttx.all_fp), 1)];

    % Raster plot
    single_unit_fig(1, example_neuron_i) = gramm('x', sdf_in.raster.stim_on([ttx.all_fp]), 'color', trial_labels);
    single_unit_fig(1, example_neuron_i).geom_raster('geom', {'point'});
    single_unit_fig(1, example_neuron_i).axe_property('XLim', xlim_vals);
    single_unit_fig(1, example_neuron_i).set_color_options('map', color_pal.foreperiod_conds(3,:));

    % Example neuron SDF
    single_unit_fig(2, example_neuron_i) = gramm('x', -2000:2000, 'y', sdf_in.sdf.stim_on([ttx.all_fp],:), 'color', trial_labels);
    single_unit_fig(2, example_neuron_i).stat_summary();
    single_unit_fig(2, example_neuron_i).axe_property('XLim', xlim_vals, 'YLim', ylim_vals{example_neuron_i});
    single_unit_fig(2, example_neuron_i).set_color_options('map', color_pal.foreperiod_conds(3,:));

end

% Add timing lines
single_unit_fig(2,1).geom_vline('xintercept',[-1300 -1000 -700 0]);
single_unit_fig(2,2).geom_vline('xintercept',[-1300 -1000 -700 0]);

% Figure setup    ////////////////////////////////////////////////////////
single_unit_fig(1,1).axe_property('XTick',[],'XColor',[1 1 1]);
single_unit_fig(1,2).axe_property('XTick',[],'XColor',[1 1 1]);

single_unit_fig(1,1).set_layout_options...
    ('Position',[0.07 0.7 0.4 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(1,2).set_layout_options...
    ('Position',[0.55 0.7 0.4 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,1).set_layout_options...
    ('Position',[0.07 0.1 0.4 0.5],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);

single_unit_fig(2,2).set_layout_options...
    ('Position',[0.55 0.1 0.4 0.5],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,...
    'margin_height',[0.00 0.00],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.0 0.00],...
    'redraw',false);


figure('Renderer', 'painters', 'Position', [100 100 600 500]);
single_unit_fig.draw();

saveas(gcf,fullfile(dirs.figures, 'example_foreperiod_all_sdf.pdf'))
