%% Tutorial for the ESCAPE EEGLAB plugin.
% 
% Copyright (C) - ESCAPE EEGLAB PLUGIN - Cedric Cannard, 2021-2025

clear; close all; clc

% launch eeglab
eeglab; close;
% pop_editoptions('option_parallel', 1); % turn parrallel computing on (1) or off (0)

% if you haven'installed the plugin yet, either go to File > Manage EEGLAB
% extensions > search ascent_compute > Install
% 
% or clone the github directory in the EEGLAB plugins folder
% or donwload the github repo and unzip it in the EEGLAB plugins folder

% Load provided sample EEG data from the tutorial directory 
% (2 minutes of resting state eyes-closed, 64-channel Biosemi):
pluginPath = fileparts(which('eegplugin_escape.m'));
cd(pluginPath)
EEG = pop_loadset('filename','sample_data.set','filepath',fullfile(pluginPath,'tutorial'));

%% Sample Entropy (SampEn) via GUI

EEG = escape_compute(EEG);  % or Tools > Compute entropy

%% SampEn via command line with default parameters

% Compute Sample entropy with command line using default parameters
EEG = escape_compute(EEG, 'measure', 'SampEn');
% print(gcf, 'SampEn_topo.png','-dpng','-r300');   % 300 dpi .png

% % Outputs can be found in:
% EEG.escape.SampEn

% Fine-tuning parameters (for demonstration only)
[EEG, com] = escape_compute(EEG, 'measure', 'SampEn', ...
    'chanlist', 'Cz Fz F3 F4 Pz Oz O1 O2', ...   % channel selection
    'tau', 2, ...
    'm', 3, ...
    'r', .5, ...
    'parallel', false, ...
    'progress', true, ...
    'vis', true);

%% Fuzzy entropy (FuzzEn)

% default parameters
EEG = escape_compute(EEG, 'measure', 'FuzzEn');

% Fine-tuning parameters (for demonstration only)
EEG = escape_compute(EEG, 'measure', 'FuzzEn', ...
    'n', 2, ...
    'kernel','gaussian', ...
    'blocksize', 128);


%% Extrema-Segmented Entropy (ExSEnt)

EEG = escape_compute(EEG, 'measure', 'ExSEnt', 'r', .15);

%% Fractal dimension (votality)

EEG = escape_compute(EEG, 'measure', 'FracDim');

%% Multiscale entropy (MSE)

EEG = escape_compute(EEG, 'measure', 'MSE', ...
    'coarsing', 'median', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 50, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'parallel', true, 'progress', true);


%% Modified MSE (mMSE)

EEG = escape_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'median', ... % 'median' (default) 'mean' 'std' 'variance'
    'num_scales', 50, ...
    'filter_mode', 'narrowband', ...  %  'narrowband' (default, annuli), 'none'
    'parallel', true, 'progress', true);

% escape_plot(EEG.escape.mMSE.data, EEG.chanlocs,'mMSE',EEG.escape.mMSE.scales)


%% Modified MSE (mMSE) - Time-resolved version

EEG = escape_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'median', ...  % 'median' (default) 'mean' 'std' 'variance'
    'num_scales', 15, ...
    'filter_mode', 'narrowband', ...  %  'narrowband' (default, annuli), 'none'
    'TimeWin', 4, ...       %  window length (in s; default = 2 for continuous data)
    'TimeStep', 2, ...    % step between centers (s); default = TimeWin/2
    'parallel', true, 'progress', true);


%% Multiscale Fuzzy Entropy (MFE)

EEG = escape_compute(EEG, 'measure', 'MFE', ...
    'coarsing', 'median', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 30, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

%% Refined Composite Multiscale Fuzzy Entropy (RCMFE)

EEG = escape_compute(EEG, 'measure', 'RCMFE', ...
    'coarsing', 'std', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 10, ...       % number of scale factors to compute (default = 20; range = 5-100 depending on sample rate)
    'n', 2, ...                 % fuzzy power (default = 2)
    'parallel', true, 'progress', true);

