function [measType, chanlist, tau, m, coarseType, nScales, filtData, n, vis, parallelComp] = ascent_compute_gui(EEG)
% Minimal wrapper that produces 3 GUIs and returns values.

measType = []; chanlist = []; tau = []; m = [];
coarseType = []; nScales = []; filtData = []; n = [];
vis = true; parallelComp = false;

% --- GUI #1: type + channels + tau + m + plot
meas = {'SampEn' 'FuzzEn' 'ExSEnt' 'FracDim' 'MSE' 'mMSE' 'MFE' 'RCMFE'};
uigeom = { [.5 .9] .5 [.5 .4 .2] .5 [.5 .1] .5 [.5 .1] .5 .5};
uilist = {
    {'style' 'text' 'string' 'Measure to compute:' 'fontweight' 'bold'}
    {'style' 'popupmenu' 'string' meas 'tag' 'etype' 'value' 5}
    {}
    {'style' 'text' 'string' 'M/EEG channels selection:' 'fontweight' 'bold'}
    {'style' 'edit' 'tag' 'chanlist'}
    {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
     'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp,tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" }
    {}
    {'style' 'text' 'string' 'Time lag (tau):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' '1' 'tag' 'tau'}
    {}
    {'style' 'text' 'string' 'Embedding dimension (m):' 'fontweight' 'bold'}
    {'style' 'edit' 'string' '2' 'tag' 'm'}
    {}
    {'style' 'checkbox' 'string' 'Plot outputs?' 'tag' 'vis' 'value' 1 'fontweight' 'bold'}
    };
param = inputgui(uigeom, uilist, 'pophelp(''ascent_compute'')','Ascent EEGLAB plugin', EEG);
if isempty(param), return; end
measType = meas{param{1}};
if ~isempty(param{2})
    chanlist = split(param{2});
else
    chanlist = {EEG.chanlocs.labels}';
end
tau  = str2double(param{3});
m    = str2double(param{4});
vis  = logical(param{5});

% --- GUI #2: coarse + nScales + filter (only for multiscale)
if contains(lower(measType), {'mse' 'mfe' 'rcmfe'})
    cTypes = {'Mean' 'Median (default)' 'Std' 'Variance'};
    uigeom = { [.5 .6] .5 [.9 .3] .5 .5 };
    uilist = {
        {'style' 'text' 'string' 'Coarse graining method:'}
        {'style' 'popupmenu' 'string' cTypes 'tag' 'stype' 'value' 2}
        {}
        {'style' 'text' 'string' 'Number of scale factors:' }
        {'style' 'edit' 'string' '20' 'tag' 'n'}
        {}
        {'style' 'checkbox' 'string' 'Filter at each scale factor (see Kosciessa et al. 2020)?','tag' 'filtscales','value',0}
        };
    param = inputgui(uigeom, uilist, 'pophelp(''ascent_compute'')','Ascent EEGLAB plugin', EEG);
    if ~isempty(param)
        coarseType = cTypes{param{1}};
        nScales    = str2double(param{2});
        filtData   = logical(param{3});
    end
end

% --- GUI #3: fuzzy power (only for fuzzy variants)
if contains(lower(measType), {'fe' 'mfe'})
    uigeom = { [.9 .3] };
    uilist = { {'style' 'text' 'string' 'Fuzzy power:' }
               {'style' 'edit' 'string' '2' 'tag' 'n'}  };
    param = inputgui(uigeom, uilist, 'pophelp(''ascent_compute'')','Ascent EEGLAB plugin', EEG);
    if ~isempty(param)
        n = str2double(param{1});
    end
end
end
