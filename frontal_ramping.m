

clear all; clc

dirs.root = '/Volumes/Mnemosyne/Codespace/2026-frontal-ramping';
dirs.raw_data = '/Volumes/Mnemosyne/Data/2026_macaque_value/spk/';

addpath(genpath(dirs.root))


%% Extract spike map
% Build (or loads if already existing) a consolidated spike_log table by 
% extracting unit-level metadata from raw .mat spike files, standardizing 
% fields, removing large waveform/timestamp variables, generating unique 
% neuron labels (session_channel), and saving the result as a CSV for 
% downstream analysis.

session_log = readtable(fullfile(dirs.root,'_data','session_log.csv'));
get_spk_map

%% Spike alignment

% Extract event-aligned sdfs and rasters
get_spk_sdfs