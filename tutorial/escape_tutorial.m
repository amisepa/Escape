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
% (several minutes of mind wandering, 64-channel Biosemi):
pluginPath = fileparts(which('eegplugin_escape.m'));
cd(pluginPath)
EEG = pop_loadset('filename','sample_data_clean.set','filepath',fullfile(pluginPath,'tutorial'));
% EEG = pop_resample(EEG, 128); % downsample to 128 Hz to increase speed
% EEG = pop_select(EEG, 'point', [1 23041]);
% EEG = ref_infinity(EEG);

%% Sample Entropy (SampEn) via GUI

EEG = escape_compute(EEG);  % or Tools > Compute entropy

%% SampEn via command line with default parameters

% Compute Sample entropy with command line using default parameters
EEG = escape_compute(EEG,'SampEn');
% print(gcf, 'SampEn_topo.png','-dpng','-r300');   % 300 dpi .png

% outputs can be found in:
EEG.escape.SampEn.data
EEG.escape.SampEn.ch

%% SampEn on just 3 channels of interest 

EEG = escape_compute(EEG,'SampEn', {'Fz' 'Cz' 'Pz'});


%% Fuzzy entropy (FuzzEn)

EEG = escape_compute(EEG,'FuzzEn');

% if you change your mind and want to see the plot, you can do so with:
escape_plot(EEG.escape.SampEn.data, EEG.escape.SampEn.electrode_locations, 'SampEn')



% Fractal votality on Fz and Cz channels
EEG = escape_compute(EEG,'Fractal votality',{'Fz' 'Cz'});

% Multiscale entropy with default parameters, on 10 time scales
EEG = escape_compute(EEG,'MSE',[],[],[],[], 10,[],[],false);
% EEG = escape_compute(EEG,'Multiscale fuzzy entropy',[],[],[],[],50,[],[],false);
escape_plot(EEG.escape.MSE.data, EEG.chanlocs, 'MFE', EEG.escape.MSE.scales)

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

