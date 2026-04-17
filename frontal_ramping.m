
%% Initialization
clear all; clc

% Define key directories
dirs.root = '/Volumes/Mnemosyne/Codespace/2026-frontal-ramping';
dirs.raw_data = '/Volumes/Mnemosyne/Data/2026_macaque_value/spk/';

% Add primary directory (and subdirectories) to path
addpath(genpath(dirs.root))

%% Setup and Information
% Load session and recording information.
session_log = readtable(fullfile(dirs.root,'_data','session_log.csv'));

% Build (or loads if already existing) a consolidated spike_log table by 
% extracting unit-level metadata from raw .mat spike files, standardizing 
% fields, removing large waveform/timestamp variables, generating unique 
% neuron labels (session_channel), and saving the result as a CSV for 
% downstream analysis.
spike_log = get_spk_map(dirs);

% Map neurons to areas and form index for future reference
area_idx = get_area_neurons(spike_log);

%% Extraction
% Extract event-aligned sdfs and rasters
% get_spk_sdfs % // ! only run to get initial extraction.


%% Ramping neuron incidence


%% Foreperiod ramping
% Separate out short, medium, long foreperiod
% (data setup: neuron x time x foreperiod condition) - averaged

foreperiod_analysis

% - Individual neurons:
%   - Condition x {slope, fano factor}

% - Population
%   - Area x foreperiod PCA


% CONSIDER FIXATION BREAKS DURING THIS PERIOD. CAN THEY BE PREDICTED FROM
% EARLIER ACTIVITY - ENGAGEMENT

%% Comparison of spiking features between areas
% ISI


%% Reward/event ramping



%% Intrinsic timescales

% https://github.com/RudebeckLab/timescales/tree/main

