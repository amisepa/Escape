%% EEGLAB plugin to compute entropy-based measures on M/EEG data.
% Also works on other types of biosignals (e.g., ECG, HRV).
% Founded on code developed by Costa et al. (2002) and Azami and Escudero
% (2017).
%
% INPUTS:
%   EEG - EEG structure in EEGLAB format
%   entropyType - 'SampEn', 'FuzzEn', 'MSE', 'MFE', 'RCMFE (default)'
%   chanlist - EEG channels of interest (empty will select all channels)
%   tau - time lag (default = 1)
%   m - embedding dimension (default = 2)
%   coarseType - coarse graining method for multiscale entropies: 'Mean',
%                 'Standard deviation' (default), or 'Variance'
%   nScales - number of scale factors (default = 30)
%   filtData - apply band pass filters to each time scale to control for
%       broadband spectral bias (see Kosciessa et al. 2020 for more detail).
%   vis - visualize entropy outputs (1, default) or not (0)
%
% USAGE:
%   1) Import EEG datra into EEGLAB
%   2) Preprocess as necessary (e.g., reference, clean data with ASR, etc.)
%   3) GUI mode: Tools > Compute entropy
% or
%   EEG = escape_compute(EEG);     % launch GUI mode
% or
%   EEG = escape_compute(EEG, 'FuzzEn',{'Fpz' 'Cz'},[],[],'Variance');
%                                       % compute fuzzy entropy only on Fpz
%                                       % and Cz channels with default tau
%                                       % and m but using variance for the
%                                       % coarse-graining
% or
%   EEG = escape_compute(EEG,'MFE',[],[],[],[],50,1,[],0);
%                                       % compute multiscale fuzzy entropy
%                                       % on all channels with default
%                                       % parameters but on 50 time scales,
%                                       % controlling for spectral bias, and
%                                       % turning plotting OFF.
%
%
% Copyright - Cedric Cannard, 2022

function [EEG, com] = escape_compute(EEG, entropyType, chanlist, tau, m, coarseType, nScales, filtData, n, vis)

com = '';

tstart =  tic;

% add path to subfolders
plugin_path = fileparts(which('eegplugin_escape.m'));
addpath(genpath(plugin_path));

% Basic checks and warnings
if nargin < 1
    help pop_entropy; return;
end
if isempty(EEG.data)
    error('Empty dataset.');
end
if isempty(EEG.chanlocs(1).labels)
    error('No channel labels.');
end
% if exist(vis,'var') && ~isempty(vis)
if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X)
    error("Electrode locations are required. " + ...
        "Go to 'Edit > Channel locations' and import the appropriate coordinates for your montage");
end
% end
if isempty(EEG.ref)
    warning(['EEG data not referenced! Referencing is highly recommended ' ...
        '(e.g., average- reference)!']);
end

% Continuous/epoched data
if length(size(EEG.data)) == 2
    continuous = true;
else
    continuous = false; %%%%%%%%%%%%% ADD OPTION TO RUN ON EPOCHED DATA %%%%%%%%%
end

%% 1st GUI to select channels and type of entropy

if nargin == 1
    eTypes = {'Sample entropy' 'Fuzzy entropy' 'Multiscale entropy' 'Multiscale fuzzy entropy' 'Refined composite multiscale fuzzy entropy'};
    uigeom = { [.5 .9] .5 [.5 .4 .2] .5 [.5 .1] .5 [.5 .1] .5 .5};
    uilist = {
        {'style' 'text' 'string' 'Entropy type:' 'fontweight' 'bold'} ...
        {'style' 'popupmenu' 'string' eTypes 'tag' 'etype' 'value' 5} ...
        {} ...
        {'style' 'text' 'string' 'Channel selection:' 'fontweight' 'bold'} ...
        {'style' 'edit' 'tag' 'chanlist'} ...
        {'style' 'pushbutton' 'string'  '...', 'enable' 'on' ...
        'callback' "tmpEEG = get(gcbf, 'userdata'); tmpchanlocs = tmpEEG.chanlocs; [tmp tmpval] = pop_chansel({tmpchanlocs.labels},'withindex','on'); set(findobj(gcbf,'tag','chanlist'),'string',tmpval); clear tmp tmpEEG tmpchanlocs tmpval" } ...
        {} ...
        {'style' 'text' 'string' 'Tau (time lag):' 'fontweight' 'bold'} ...
        {'style' 'edit' 'string' '1' 'tag' 'tau'} ...
        {} ...
        {'style' 'text' 'string' 'Embedding dimension:' 'fontweight' 'bold'} ...
        {'style' 'edit' 'string' '2' 'tag' 'm'}  ...
        {} ...
        {'style' 'checkbox' 'string' 'Plot entropy outputs' 'tag' 'vis' 'value' 1 'fontweight' 'bold'}  ...
        };
    param = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
    entropyType = eTypes{param{1}};
    if ~isempty(param{2})
        chanlist = split(param{2});
    else
        chanlist = {EEG.chanlocs.labels}';
    end
    tau  = str2double(param{3});
    m = str2double(param{4});
    vis = logical(param{5});
end


%% 2nd GUI to select additional parameters

if nargin == 1 && contains(lower(entropyType), 'multiscale')
    cTypes = {'Mean' 'Standard deviation (default)' 'Variance'};
    uigeom = { [.5 .6] .5 [.9 .3] .5 .5 };
    uilist = {
        {'style' 'text' 'string' 'Coarse graining method:'} ...
        {'style' 'popupmenu' 'string' cTypes 'tag' 'stype' 'value' 2} ...
        {} ...
        {'style' 'text' 'string' 'Number of scale factors:' } ...
        {'style' 'edit' 'string' '30' 'tag' 'n'}  ...
        {} ...
        {'style' 'checkbox' 'string' 'Bandpass filter each time scale (recommended to control for spectral bias)','tag' 'filtData','value',0}  ...
        };
    param = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','get_entropy() EEGLAB plugin',EEG);
    coarseType = cTypes{param{1}};
    nScales = str2double(param{2});
    filtData = logical(param{3});
end

%% 3rd GUI for fuzzy power

if nargin == 1 && contains(entropyType, 'Fuzzy')
    uigeom = { [.9 .3] };
    uilist = { {'style' 'text' 'string' 'Fuzzy power:' } ...
        {'style' 'edit' 'string' '2' 'tag' 'n'}  };
    param = inputgui(uigeom,uilist,'pophelp(''pop_entropy'')','entropy EEGLAB plugin',EEG);
    n = str2double(param{1});
end

%% Defaults if something was missed in command line

if ~exist('chanlist','var') || isempty(chanlist)
    disp('No channels were selected: selecting all channels (default)')
    chanlist = {EEG.chanlocs.labels}';
end
if ~exist('entropyType','var') || isempty(entropyType)
    disp('No entropy type selected: selecting Refined composite multiscale fuzzy entropy (default)')
    entropyType = 'Refined composite multiscale fuzzy entropy';
end
if ~exist('tau','var') || isempty(tau)
    disp('No time lag selected: selecting tau = 1 (default).')
    tau = 1;
end
if ~exist('m','var') || isempty(m)
    disp('No embedding dimension selected: selecting m = 2 (default).')
    m = 2;
end
if ~exist('vis','var') || isempty(vis)
    disp('Plotting option not selected: turning plotting ON (default).')
    vis = true;
end
if contains(lower(entropyType), {'mse' 'mfe' 'rcmfe'})
    if ~exist('coarseType','var') || isempty(coarseType)
        disp('No coarse graining method selected: selecting standard deviation (default).')
        coarseType = 'Standard deviation';
    end
    if ~exist('nScales','var') || isempty(nScales)
        disp('Number of scale factors not selected: selecting nScales = 30 (default).')
        nScales = 30;
    end
    if ~exist('filtData','var') || isempty(filtData)
        filtData = false;
    end
else
    coarseType = [];
    nScales = [];
    filtData = [];
end
if contains(lower(entropyType), {'fuzzen' 'mfe' 'rcmfe'})
    if ~exist('n','var') || isempty(n)
        disp('No fuzzy power selected: selecting n = 2 (default).')
        n = 2;
    end
else
    n = [];
end

% r = .15; % Hardcode r to .15 because data are later normalized to have SD of 1 (Azami approach)

% % parallel computing
parallelComp = false;   % request parallel mode (to add to GUI and command line options)
p = gcp('nocreate');   % get current pool (if any)
if parallelComp
    if isempty(p)
        n = max(2, min(6, feature('numcores')-1));  % safe
        % n = feature('numcores')-1; % use almost all cores
        p = parpool('Processes', n);
        pctRunOnAll maxNumCompThreads(1);   % keep each worker single-threaded
        % pctRunOnAll setenv('OMP_NUM_THREADS','1'); setenv('MKL_NUM_THREADS','1');
        fprintf('Started new pool with %d workers.\n', p.NumWorkers);
    else
        fprintf('Pool already running with %d workers.\n', p.NumWorkers);
    end
else
    if ~isempty(p)
        delete(p);
        fprintf('Closed existing pool.\n');
    end
end

%% Compute entropy depending on choices

% index with channels of interest
nChan = length(chanlist);
if nChan > 1 && nChan < EEG.nbchan
    [~, chanIdx] = intersect({EEG.chanlocs.labels}, split(chanlist));
else
    chanIdx = 1:EEG.nbchan;
end

% extract data of interest and avoid structure broadcast overhead for
% parfor loops
data = EEG.data(chanIdx, :);    % EEG data (channels x samples)
fs = EEG.srate;                 % sample rate
nChan = size(data, 1);          % number of channels
chanLabels = {EEG.chanlocs(chanIdx).labels}; % channel labels
chanlocs = EEG.chanlocs(chanIdx);

% preallocate memory
if contains(lower(entropyType), {'mse' 'mfe' 'rcmfe'})
    entropy = nan(nChan, nScales);
else
    entropy = nan(nChan,1);
end
scales = {};

% COMPUTE ENTROPY/COMPLEXITY MEASURES
switch entropyType

    % Sample Entropy (SampEn)
    case 'SampEn'
        tic
        entropy = compute_SampEn(data, 'm', m, 'tau', tau, 'Parallel', parallelComp);
        toc
        
    % Fuzzy Entropy (FuzzEn)
    case 'FuzzEn'
        
        % Compute fuzzy entropy per channel
        kernel_meth = 'exponential'; % 'exponential' (default) or 'gaussian'
        entropy = compute_FuzzEn(data, 'm', m, 'tau', tau, 'n', n, ...
            'Kernel', kernel_meth, 'Parallel', parallelComp);
        toc


    % Fractal Dimension (FracDim)
    case 'FracDim'
        entropy = fractal_volatility(data);

    case 'MSE'
        % [entropy, scales, info] = compute_mse(data, fs, m, r, tau, coarseType, nScales, filtData);

        %   'FilterMode'       : 'lowpass' (default) | 'highpass' | 'bandpass' | 'bandstop' | 'none'
        [entropy, scales, info] = compute_mse(data, fs, ...
            m, r, tau, coarseType, nScales, filtData, ...
            'FilterMode','bandpass', 'BandMode','annuli', ...   % <- for non-overlapping rings
            'CoarseMethod','skip', ...                          % <- true filt-skip (optional)
            'Parallel', parallelComp, 'Progress',true, 'ProgressStyle','waitbar');

        % % Classic mode
        % [entropy, scales, info] = compute_mse(data, fs, ...
        %     m, r, tau, coarseType, nScales, true, ...
        %     'FilterMode','lowpass', 'CoarseMethod','stat', ...
        %     'Parallel','none', 'Verbose',true, 'Progress',true, 'ProgressStyle','waitbar');


    case 'MFE'
        % disp('Computing multiscale fuzzy entropy...')
        % progressbar('Channels')
        % % t1 = tic;
        % for iChan = 1:nChan
        %     fprintf('Channel %d: \n', iChan)
            [entropy(iChan,:), scales] = compute_mfe(data, ...
                m, r, tau, coarseType, nScales, filtData, EEG.srate, n);
            % progressbar(iChan/nChan)
            % % entropy(iChan,1:length(enttmp)) = enttmp;
        % end
        % toc(t1)

        % % Remove NaN scales
        % idx = isnan(entropy(1,:));
        % entropy(:,idx) = [];
        % scales(idx) = [];

    case 'RCMFE'
        % disp('Computing multiscale fuzzy entropy...')
        % progressbar('Channels')
        % for iChan = 1:nChan
            % fprintf('Channel %d: \n', iChan)
            % signal = EEG.data(chanIdx(iChan),:);
            [entropy(iChan,:), scales] = compute_rcmfe(signal, ...
                m, r, tau, coarseType, nScales, filtData, EEG.srate, n);
            % [entropy(iChan,:), scales] = compute_rcmfe(EEG.data(iChan,:),m,r,tau,coarseType,nScales,filtData,EEG.srate);
            % progressbar(iChan/nChan)
        % end

    % Spectral entropy (SpecEn) - Over time
    case 'SpecEn' 
        disp('Computing spectral entropy...')
        progressbar('Channels')
        parfor iChan = 1:nchan
            % fprintf('Channel %d \n', iChan)
            [entropy(iChan,:),te] = pentropy( zscore(EEG.data(chanIdx(iChan),:)), EEG.srate);

            fprintf('   %s: %6.3f \n', EEG.chanlocs(iChan).labels, entropy(iChan,:))
            progressbar(iChan/nchan)
        end

    otherwise
        error('Unknown entropy type. Please select one of the options (see help get_entropy).')
end

% % Get scales bounds
% if ~isempty(scales)
%     scale_bounds = {};
%     for iScale = 1:length(scales)
%         % Scale frequency bounds
%         nf = EEG.srate / 2; % Nyquist frequency
%         upperBound = (1/iScale).*nf + .05*((1./iScale).*nf);
%         lowerBound = (1/(iScale+1)).*nf - .05*((1./(iScale+1)).*nf);
%         scale_bounds{iScale} = [round(upperBound,3) round(lowerBound,3) ];
%     end
% end

% Visualize
if vis
    if nChan>1
        escape_plot(entropy, chanlocs, entropyType, scales);
    else
        disp("Sorry, you need more than 1 EEG channel for visualization.")
    end
end

% save outputs in EEG structure
EEG.escape.(entropyType).data = entropy;
EEG.escape.(entropyType).electrode_labels = chanLabels;
EEG.escape.(entropyType).electrode_locations = chanlocs;
if contains(lower(entropyType), {'mse' 'mfe' 'rcmfe'})
    EEG.escape.(entropyType).scales = scales;
end

%%%%%% ADD PARAMETERS USED IN STRUCTURE OUTPUT %%%%%%%%%%


% Command history (TO FIX)
chanLabels = strjoin(chanlist);
chanLabels = insertBefore(chanLabels," ", "'");
chanLabels = insertAfter(chanLabels," ", "'");
com = sprintf('EEG = escape_compute(''%s'', {''%s''}, %d, %d, %s, %d, %d, %s, %d);', ...
    entropyType,chanLabels,tau,m,coarseType,nScales,filtData,'[]',vis);

gong
disp('Done computing with Escape! Outputs can be found in the EEG.escape structure.')
fprintf('Time to compute: %.2f minutes. \n', toc(tstart)/60)

