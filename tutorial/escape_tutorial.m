%% Tutorial for the ASCENT EEGLAB plugin.
% 
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2021-2025

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
    'coarsing', 'std', ...     % 'median' (default) 'mean' 'trimmed mean' 'std' 'var'
    'num_scales', 30, ...       % number of scale factors to compute (default = 15; range = 5-100 depending on sample rate)
    'parallel', true, 'progress', true);

%% Modified MSE (mMSE)

EEG = escape_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'median', ... % 'median' (default) 'mean' 'std' 'variance'
    'num_scales', 20, ...
    'filter_mode', 'narrowband', ...  %  'narrowband' (default, annuli), 'none'
    'parallel', true, 'progress', true);

%% Modified MSE (mMSE) - Time-resolved version

EEG = escape_compute(EEG, 'measure', 'mMSE', ...
    'coarsing', 'std', ...  % 'median' (default) 'mean' 'std' 'variance'
    'num_scales', 20, ...
    'filter_mode', 'narrowband', ...  %  'narrowband' (default, annuli), 'none'
    'TimeWin', 2, ...       %  window length (in s; default = 2 for continuous data)
    'TimeStep', 0.5, ...    % step between centers (s); default = TimeWin/2
    'parallel', true, 'progress', true);


%% Multiscale fuzzy entropy (MFE)

coarsing = 'median';  % 'median' 'mean' 'sd' 'variance'
EEG = escape_compute(EEG,'MFE',[],[],[], coarsing, 10, false);


%% Multiscale fuzzy entropy (MFE)

coarsing = 'median';  % 'median' 'mean' 'sd' 'variance'
EEG = escape_compute(EEG,'RCMFE',[],[],[], coarsing, 10, false);




%% 
% % if you change your mind and want to see the plot, you can do so with:
% escape_plot(EEG.escape.SampEn.data, EEG.escape.SampEn.electrode_locations, 'SampEn')

% same but only on channels F3 and F4 (and plotting On)
EEG = escape_compute(EEG,'MFE',{'F3' 'F4'},[],[],[],10,[],[],true,false);

% Same but using bandpass filtering at each time scale to control for
% spectral contamination (and plotting On)
EEG = escape_compute(EEG,'MFE',{'F3' 'F4'},[],[],[],10,true,[],true);

% EEG = escape_compute(EEG);   % GUI mode
% EEG = escape_compute(EEG,'Approximate entropy');
% EEG = escape_compute(EEG,'Sample entropy');
% EEG = escape_compute(EEG,'Fuzzy entropy',{'Cz'});
% EEG = escape_compute(EEG,'Multiscale entropy', {'Cz' 'O1' 'Fz'}, [],[],[],30);
% EEG = escape_compute(EEG,'Multiscale fuzzy entropy', {'Cz' 'O1' 'Fz'}, [],[],[],30);
% EEG = escape_compute(EEG, 'Sample entropy',[], 1, 2, 'Mean', 1, 2);
% EEG = escape_compute(EEG, 'Multiscale entropy', {EEG.chanlocs.labels}, 1, 2, 'Standard deviation',15,1,[],1);

