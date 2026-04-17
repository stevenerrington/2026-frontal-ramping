function spike_log = get_spk_map(dirs)

% Construct output file path once
outFile = fullfile(dirs.root,'_data','spike_log.csv');


if exist(outFile, 'file') == 2
    % ---- IMPORT PRECREATED CSV ----
    spike_log = readtable(outFile);
else
    % ---- BUILD FROM RAW FILES ----
    % Get list of .mat files
    files = dir('/Volumes/Mnemosyne/Data/2026_macaque_value/spk/*.mat');
    fileNames = {files.name};

    % Preallocate container for tables (faster than growing table each loop)
    all_rows = cell(numel(fileNames), 1);
    row_counter = 0;

    % Loop over files
    for file_i = 1:numel(fileNames)

        % Load file
        data_in = load(fullfile(dirs.raw_data, fileNames{file_i}));

        % Skip if 'unit' field is missing or empty
        if ~isfield(data_in, 'unit') || isempty(data_in.unit)
            continue
        end

        units = data_in.unit;
        unit_tables = cell(numel(units), 1);

        % Loop over units (neurons)
        for unit_i = 1:numel(units)
            s = units{unit_i};

            % ---- CLEAN FIELDS ----

            % Ensure 'area' is string
            if isfield(s, 'area') && ischar(s.area)
                s.area = string(s.area);
            end

            % Store timestamps as cell (for table compatibility)
            if isfield(s, 'timestamps')
                s.timestamps = {s.timestamps};
            end

            % Convert struct to table row
            unit_tables{unit_i} = struct2table(s);
        end

        % Concatenate all units from this file once
        file_table = vertcat(unit_tables{:});

        % Store in master list
        row_counter = row_counter + 1;
        all_rows{row_counter} = file_table;
    end

    % Remove unused preallocated cells
    all_rows = all_rows(1:row_counter);

    % Concatenate everything once (FAST)
    spike_log = vertcat(all_rows{:});

    % ---- POST-PROCESSING ----
    % Remove unwanted variables (do it in one call)
    vars_to_remove = intersect( ...
        ["avg_waveform","std_waveform","timestamps"], ...
        spike_log.Properties.VariableNames);

    spike_log = removevars(spike_log, vars_to_remove);

    % Create neuron label: session + "_" + channel
    spike_log.neuron_label = spike_log.session + "_" + spike_log.ch + "_" + spike_log.clust_id;

    % ---- SAVE ----
    writetable(spike_log, outFile, 'WriteRowNames', true);

end
end