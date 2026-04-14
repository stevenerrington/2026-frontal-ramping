

%% SDF & Raster Extraction -----------------------------------------------------------------
% Set convolution parameters
ops.timewin = -2000:2000;
ops.sdf_filter = 'PSP';

% Loop through sessions
parfor session_i = 1:size(session_log,1)
    disp(['Loading session ' int2str(session_i) ' of ' int2str(size(session_log,1))])

    data_in = load(fullfile(dirs.raw_data,[session_log.session{session_i}, '_spk.mat']));

    % Align spikes and generate SDF/rasters
    for neuron_i = 1:size(data_in.unit,2)
        disp(['   - processing neuron ' int2str(neuron_i) ' of ' int2str(size(data_in.unit,2))])
        spikes = [];
        spikes = data_in.unit{neuron_i}.timestamps * 1000; % Get unit spike times and convert to milliseconds

        sdf = struct(); raster = struct();

        [sdf.fixation, raster.fixation] = get_aligned_sdf(spikes, round(data_in.t_evt.fixcross_fix)*1000, ops); % Align units to trial event (fixation)
        [sdf.stim_on, raster.stim_on] = get_aligned_sdf(spikes, round(data_in.t_evt.stim_on)*1000, ops); % Align units to trial event (target)
        [sdf.reward, raster.reward] = get_aligned_sdf(spikes, round(data_in.t_evt.rew_on)*1000, ops); % Align units to trial event (reward)
        [sdf.brk_fix, raster.brk_fix] = get_aligned_sdf(spikes, round(data_in.t_evt.brk_fix)*1000, ops); % Align units to trial event (break fixation)


        data_out = struct('sdf', sdf, 'raster', raster);
        out_filename = [data_in.unit{neuron_i}.session '_' data_in.unit{neuron_i}.ch '.mat'];
        
        save(fullfile(dirs.root,'_data','sdf',out_filename),'-fromstruct', data_out, '-v7.3')

    end
end