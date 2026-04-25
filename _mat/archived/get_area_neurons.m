function area_idx = get_area_neurons(spike_log)
% Find all areas within dataset
areas = unique(spike_log.area);

% Loop through areas and map neuron idx to recording area
area_idx = struct();

for area_i = 1:length(areas)
    area_idx.(['a_' areas{area_i}]) = find(strcmp(spike_log.area, areas{area_i}));
end
