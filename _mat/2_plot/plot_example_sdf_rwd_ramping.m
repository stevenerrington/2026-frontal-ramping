%% Example foreperiod-ramping neuron

% Select an example neuron
example_neuron_idx = neuron_idx.rwd_pos_ramping(21);

% Load neuron and session data
neuron_label = spike_log.neuron_label{example_neuron_idx};

beh_data = load(fullfile( ...
    dirs.raw_data, ...
    [spike_log.session{example_neuron_idx}, '_spk.mat']));

sdf_data = load(fullfile( ...
    dirs.proc_data, ...
    [neuron_label, '.mat']));

%% Select trials reward

% Restrict to trials with approximately 1-second foreperiod
trial_mask = beh_data.trialtype.rew == 1;

% Extract rasters and spike-density functions aligned to stimulus onset
example_raster = sdf_data.raster.stim_on(trial_mask);
example_sdf    = sdf_data.sdf.reward(trial_mask, :);

%% Smooth single-trial SDFs

for trial = 1:size(example_sdf, 1)
    example_sdf(trial, :) = smooth(example_sdf(trial, :), 100)';
end

%% Plot settings

time_axis = -2000:2000;      % ms relative to stimulus onset
x_limits  = [-1250 250];

clear single_unit_fig

%% Raster plot

single_unit_fig(1,1) = gramm('x', example_raster);

single_unit_fig(1,1).geom_raster('geom', {'point'});
single_unit_fig(1,1).axe_property('XLim', x_limits);

%% Population SDF

single_unit_fig(2,1) = gramm( ...
    'x', time_axis, ...
    'y', example_sdf);

single_unit_fig(2,1).stat_summary();
single_unit_fig(2,1).axe_property( ...
    'XLim', x_limits, ...
    'YLim', [15 45]);

%% Layout configuration

% Raster panel
single_unit_fig(1,1).set_layout_options( ...
    'Position', [0.1 0.75 0.8 0.2], ...
    'legend', false, ...
    'margin_height', [0.00 0.00], ...
    'margin_width', [0.00 0.00], ...
    'redraw', false);

% SDF panel
single_unit_fig(2,1).set_layout_options( ...
    'Position', [0.1 0.2 0.8 0.5], ...
    'legend', false, ...
    'margin_height', [0.00 0.00], ...
    'margin_width', [0.00 0.00], ...
    'redraw', false);

%% Axis formatting and event markers

% Hide x-axis from raster panel
single_unit_fig(1,1).axe_property( ...
    'XTick', [], ...
    'XColor', [1 1 1]);

% Mark foreperiod start (-1000 ms) and stimulus onset (0 ms)
single_unit_fig(2,1).geom_vline( ...
    'xintercept', [-1000 0], ...
    'style', {'k-','k-'});

%% Draw figure

figure( ...
    'Position', [100 100 300 420], ...
    'Renderer', 'painters');

single_unit_fig.draw();