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
pluginPath = fileparts(which('eegplugin_Ascent.m'));
cd(pluginPath)
EEG = pop_loadset('filename','sample_data_clean.set','filepath',fullfile(pluginPath,'tutorial'));
% EEG = pop_resample(EEG, 128); % downsample to 128 Hz to increase speed
% EEG = pop_select(EEG, 'point', [1 23041]);
% EEG = ref_infinity(EEG);

% Launch GUI to selec all parameters manually
EEG = ascent_compute(EEG);  % or Tools > Compute entropy

% Compute Sample entropy with command line using default parameters
EEG = ascent_compute(EEG,'Sample entropy', {EEG.chanlocs(1:3:end).labels});
print(gcf, 'SampEn_topo.png','-dpng','-r300');   % 300 dpi .png

% If you forgot to plot outputs and want to see, you can call the function
% like this: 
ascent_plot(EEG.entropy, EEG.chanlocs, 1:EEG.nbchan)

% same but only on Fz and Cz channels and using variance instead of default sd
EEG = ascent_compute(EEG,'Fuzzy entropy',{EEG.chanlocs(1:10:end).labels},[],[],'Variance');

% Fractal votality on Fz and Cz channels
EEG = ascent_compute(EEG,'Fractal votality',{'Fz' 'Cz'});

% Multiscale entropy with default parameters, on 10 time scales
EEG = ascent_compute(EEG,'Multiscale entropy',[],[],[],[],10,[],[],false);
% EEG = ascent_compute(EEG,'Multiscale fuzzy entropy',[],[],[],[],50,[],[],false);

% same but only on channels F3 and F4 (and plotting On)
EEG = ascent_compute(EEG,'Multiscale fuzzy entropy',{'F3' 'F4'},[],[],[],10,[],[],true,false);

% Same but using bandpass filtering at each time scale to control for
% spectral contamination (and plotting On)
EEG = ascent_compute(EEG,'Multiscale fuzzy entropy',{'F3' 'F4'},[],[],[],10,true,[],true);

% EEG = ascent_compute(EEG);   % GUI mode
% EEG = ascent_compute(EEG,'Approximate entropy');
% EEG = ascent_compute(EEG,'Sample entropy');
% EEG = ascent_compute(EEG,'Fuzzy entropy',{'Cz'});
% EEG = ascent_compute(EEG,'Multiscale entropy', {'Cz' 'O1' 'Fz'}, [],[],[],30);
% EEG = ascent_compute(EEG,'Multiscale fuzzy entropy', {'Cz' 'O1' 'Fz'}, [],[],[],30);
% EEG = ascent_compute(EEG, 'Sample entropy',[], 1, 2, 'Mean', 1, 2);
% EEG = ascent_compute(EEG, 'Multiscale entropy', {EEG.chanlocs.labels}, 1, 2, 'Standard deviation',15,1,[],1);

