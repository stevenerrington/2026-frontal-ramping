
%% Initialization
clear all; clc

% Define key directories
dirs.root = '/Volumes/Mnemosyne/Codespace/2026-frontal-ramping';
dirs.raw_data = '/Volumes/Mnemosyne/Data/2026_macaque_value/spk/';
dirs.proc_data = '/Volumes/Mnemosyne/Data/2026_macaque_value/proc/';
dirs.figures = '/Volumes/Mnemosyne/Codespace/2026-frontal-ramping/_figures/';

% Add primary directory (and subdirectories) to path
addpath(genpath(dirs.root))

%% ========================================================================
% PREPROCESSING (MATLAB)
%
% MATLAB is now used exclusively to extract the raw arrays that every
% downstream analysis is built on. All analysis, statistics, and
% visualisation happens in Python — see frontal_ramping_analysis.ipynb.
% =========================================================================

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

%% Extraction: event-aligned SDFs and rasters
% Adapted to run on PC cluster (see _mat/adapt/get_spk_sdfs_PCADAPTER.m).
% Data is then transferred back locally for the array extraction step
% below. Also extracts the per-neuron average waveform array.

%% Extraction: condition-split trial-averaged arrays
% Loops through every neuron once and produces the three arrays that
% feed all downstream Python analysis:
%   fixation_array.mat — stim-onset-aligned SDF, split by foreperiod
%   fixon_array.mat    — fixation-onset-aligned SDF, split by foreperiod
%   reward_array.mat   — reward-aligned SDF, split by chosen reward probability
extract_arrays

%% ========================================================================
% ANALYSIS & VISUALISATION (Python)
%
% Everything downstream of the three extracted arrays — z-scoring, ramp
% detection, prevalence/odds-ratio statistics, SDF and PCA trajectory
% figures, ramp-onset latency estimation — lives in
% frontal_ramping_analysis.ipynb.
% =========================================================================
