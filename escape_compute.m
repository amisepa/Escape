function [EEG, com] = escape_compute(EEG, varargin)

%% EEGLAB plugin to compute entropy-based measures on M/EEG data.
% Also works on other types of biosignals (e.g., ECG, HRV).
% Founded on code developed by Costa et al. (2002) and Azami and Escudero
% (2017).
%
% INPUTS:
%   EEG - EEG structure in EEGLAB format
%   measure - 'SampEn', 'FuzzEn', 'MSE', 'MFE', 'RCMFE (default)'
%   chanlist - EEG channels of interest (empty will select all channels)
%   tau - time lag (default = 1)
%   m - embedding dimension (default = 2)
%   coarsing - coarse graining method for multiscale entropies: 'Mean',
%                 'Standard deviation' (default), or 'Variance'
%   num_scales - number of scale factors (default = 30)
%   filt_scales - apply band pass filters to each time scale to control for
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
%   EEG = escape_compute(EEG,'MFE',[],[],[],[],50,true,[],0);
%                                       % compute multiscale fuzzy entropy
%                                       % on all channels with default
%                                       % parameters but on 50 time scales,
%                                       % controlling for spectral bias, and
%                                       % turning plotting OFF.
%
%
% Copyright - Cedric Cannard, 2022


tstart = tic;

% add path to subfolders
plugin_path = fileparts(which('eegplugin_escape.m'));
addpath(genpath(plugin_path));

% Basic checks and warnings (unchanged)
if nargin < 1, help pop_entropy; return; end
if isempty(EEG.data), error('Empty dataset.'); end
if isempty(EEG.chanlocs(1).labels), error('No channel labels.'); end
if ~isfield(EEG.chanlocs, 'X') || isempty(EEG.chanlocs(1).X)
    error("Electrode locations are required. Go to 'Edit > Channel locations' and import coordinates");
end
if isempty(EEG.ref)
    warning('EEG data not referenced! Referencing is highly recommended (e.g., average-reference)!');
end

% Continuous/epoched data (unchanged)
if length(size(EEG.data)) == 2
    continuous = true;
else
    continuous = false; %%%%%%%%%%%%% ADD OPTION TO RUN ON EPOCHED DATA %%%%%%%%%
end

measure     = [];
chanlist    = [];
tau         = [];
m           = [];
coarsing    = [];
num_scales  = [];
filt_scales = [];
n           = [];
vis         = true;
paraComp    = false;

% -----------------------------
% Get params: GUI if none, else name–value pairs
% -----------------------------
if nargin == 1 || (nargin == 2 && isempty(varargin{1}))
    [measure, chanlist, tau, m, coarsing, num_scales, filt_scales, n, vis, paraComp] = ...
        escape_compute_gui(EEG);
else
    if mod(numel(varargin),2) ~= 0
        error('Options must be provided as name–value pairs.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k+1};
        switch key
            case 'measure'
                measure = val;
            case 'chanlist'
                if ischar(val) || isstring(val)
                    % allow "all" or space/comma-separated labels (same behavior as GUI)
                    if strcmpi(string(val),'all')
                        chanlist = {EEG.chanlocs.labels}';
                    else
                        chanlist = regexp(char(val),'[^,\s]+','match')';
                    end
                elseif iscell(val)
                    chanlist = val(:);
                else
                    error('Channels must be ''all'', a string list, or a cellstr.');
                end
            case 'tau'
                tau = double(val);
            case 'm'
                m = double(val);
            case 'r'
                r = double(val);
            case 'coarsing'
                coarsing = val;
            case 'num_scales'
                num_scales = double(val);
            case 'filt'
                filt_scales = logical(val);
            case 'n'
                n = double(val);
            case 'vis'
                vis = logical(val);
            case 'parallel'
                paraComp = logical(val);
            case 'progress'
                trackProg = logical(val);
            case 'kernel'
                kernel_meth = val;
            case 'blocksize'
                blocksize = double(val);
            otherwise
                error('Unknown option: %s', key);
        end
    end
end

% --------
% Defaults
% --------

% general parameters
if ~exist('chanlist','var') || isempty(chanlist)
    disp('No channels were selected: selecting all channels (default)')
    chanlist = {EEG.chanlocs.labels}';
end
if ~exist('measure','var') || isempty(measure)
    disp('No entropy type selected: selecting Refined composite multiscale fuzzy entropy (default)')
    measure = 'Refined composite multiscale fuzzy entropy';
end
if ~exist('tau','var') || isempty(tau)
    disp('No time lag selected: selecting tau = 1 (default).')
    tau = 1;
end
if ~exist('m','var') || isempty(m)
    disp('No embedding dimension selected: selecting m = 2 (default).')
    m = 2;
end
if ~exist('r','var') || isempty(r)
    disp('No similarity bound selected: selecting r = .15 (default).')
    r = .15;
end
if ~exist('vis','var') || isempty(vis)
    disp('No plotting option not selected: turning plotting ON (default).')
    vis = true;
end
if ~exist('progress','var') || isempty(trackProg)
    disp('No progress tracking defined: setting it to ON (default).')
    trackProg = true;
end
if ~exist('parallel','var') || isempty(paraComp)
    disp('Computing method not selected: turning parallel computing ON (default).')
    paraComp = true;
end

% multiscale parameters
if contains(lower(measure), {'mse' 'mfe' 'rcmfe'})
    if ~exist('coarsing','var') || isempty(coarsing)
        disp('No coarse graining method selected: selecting standard deviation (default).')
        coarsing = 'Standard deviation';
    end
    if ~exist('num_scales','var') || isempty(num_scales)
        disp('Number of scale factors not selected: selecting num_scales = 30 (default).')
        num_scales = 30;
    end
    if ~exist('filt_scales','var') || isempty(filt_scales)
        filt_scales = false;
        disp("Filtering each scale factor (from Kosciessa et al. 2017) set to: OFF.")
    else
        disp("Filtering each scale factor (from Kosciessa et al. 2017) set to: ON.")
    end
else
    coarsing = [];
    num_scales = [];
    filt_scales = [];
end

% Fuzzy parameters
if contains(lower(measure), {'fuzzen' 'mfe' 'rcmfe'})
    if ~exist('n','var') || isempty(n)
        disp('No fuzzy power selected: using n = 2 (default).')
        n = 2;
    end
    if ~exist('kernel','var') || isempty(kernel_meth)
        disp('Kernel method not selected: using exponential kernel (default).')
        kernel_meth = 'exponential';
    end
    if ~exist('blocksize','var') || isempty(blocksize)
        disp('BlockSize not specified: using blocksize = 256 (default).')
        blocksize = 256;
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
if contains(lower(measure), {'mse' 'mfe' 'rcmfe'})
    entropy = nan(nChan, num_scales);
else
    entropy = nan(nChan,1);
end
scales = {};

% COMPUTE ENTROPY/COMPLEXITY MEASURES
switch measure

    % Sample Entropy (SampEn)
    case 'SampEn'
        % entropy = compute_SampEn(data, 'm', m, 'tau', tau, 'Parallel', parallelComp);
        entropy = compute_SampEn(data, 'm', m, 'r', r, 'parallel', paraComp, 'progress', trackProg);

        
    % Fuzzy Entropy (FuzzEn)
    case 'FuzzEn'
        % kernel_meth = 'exponential'; % 'exponential' (default) or 'gaussian'
        entropy = compute_FuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
            'Kernel', kernel_meth, 'BlockSize',blocksize, 'Parallel', paraComp, 'progress', trackProg);

    % Extrema-Segmented Entropy (ExSEnt)    
    case 'ExSEnt'
        % HD = nan(nChan,1);
        % HA = nan(nChan,1);
        % HDA = nan(nChan,1);
        % parfor iChan = 1:nChan
        %     fprintf('channel %g/%g\n', iChan, nChan)
        %     [HD(iChan,:), HA(iChan,:), HDA(iChan,:)] = compute_ExSEnt(data(iChan,:), 0.001, m, r);
        % end

        [HD, HA, HDA, info] = compute_ExSEnt2(data, 'm', m, ...
            'r', r, 'lambda', 0.001, 'Plot', false, ...  % add these options to GUI/command?
            'Parallel', paraComp);

    % Fractal Dimension (FracDim)
    case 'FracDim'
        progressbar('Computing fractal dimension on all channels')
        % for iChan = 1:nChan
        %     [entropy(iChan,:), sd] = fractal_volatility(data(iChan,:)');
        %     progressbar(iChan / nChan);
        % end
        [entropy, SD, info] = compute_FracDim(data, ...
            'RejectBursts', true, 'WinFrac', 0.02, 'ZThresh', 6, ...
            'RobustFit', 'ols', ... %  'theilsen' | 'ols'
            'ScaleTrimIQR', true, 'MinJ', [], 'MaxJ', [], ...
            'Parallel', paraComp);

    case 'MSE'

        % % Origianl Costa 2002 MSE
        % parfor iChan = 1:nChan
        %     fprintf('channel %g/%g\n',iChan, nChan);
        %     % [ entropy(iChan,:), A, B ] = mse(data(iChan,:));
        %     entropy(iChan,:) = compute_mse_costa(data(iChan,:), m, .15, tau, 'mean', num_scales, filt_scales, fs)
        % end
        % if ~filt_scales
        %     scales = arrayfun(@(x) {num2str(x)}, 1:num_scales);
        % end
        
      % Enhanced classic MSE
      % [entropy, scales] = compute_mse(EEG, 'Fs', fs, 'm', m, 'tau', tau, ...
      %     'FilterMode', 'none', 'coarsing', coarsing, 'num_scales', num_scales, ...
      %     'Parallel', parallelComp, 'Progress', true);           
         
         if filt_scales 
             % Kosciessa-style narrowband MSE (default): annuli passbands + filter-skip
             filter_mode = 'narrowband';  %  'narrowband' (default, annuli) | 'lowpass' | 'highpass' | 'none'
         else
             filter_mode = 'none';
         end

         [entropy, scales, info] = compute_mse(EEG, 'Fs', fs, 'm', m, 'tau', tau, 'r', r, ...
              'FilterMode', filter_mode, 'coarsing', coarsing, 'num_scales', num_scales, ...
              'Parallel', paraComp, 'Progress', true);

          % % Time-resolved version
          % [entropy, scales, info] = compute_mse(EEG, 'Fs', fs, 'm', m, 'tau', tau, ...
          %      'FilterMode', filter_mode, 'coarsing', coarsing, 'num_scales', num_scales, ...
          %     'TimeWin', 2, 'TimeStep', 1, ...
          %     'Parallel', parallelComp, 'Progress', true);

          
    case 'MFE'
        % [entropy, scales] = compute_mfe(data, ...
        %     m, r, tau, coarsing, num_scales, filt_scales, EEG.srate, n);

         filter_mode = 'narrowband';  %  'narrowband' (default, annuli) or 'none'
         [entropy, scales, info] = compute_mfe(EEG, 'Fs', fs, 'm', m, 'tau', tau, 'r', r, ...
              'FilterMode', filter_mode, 'coarsing', coarsing, 'num_scales', num_scales, ...
              'Parallel', paraComp, 'Progress', true);

          % % Time-resovled version
          % [entropy, scales, info] = compute_mse(EEG, 'Fs', fs, 'm', m, 'tau', tau, ...
          %      'FilterMode', filter_mode, 'coarsing', coarsing, 'num_scales', num_scales, ...
          %     'TimeWin', 2, 'TimeStep', 1, ...
          %     'Parallel', parallelComp, 'Progress', true);


    case 'RCMFE'
        parfor iChan = 1:nChan
            fprintf('channel %g/%g\n',iChan, nChan);
            entropy(iChan,:) = compute_rcmfe(data(iChan,:), ...
                m, .15, tau, coarsing, num_scales, fs, n, false);
        end
        scales = arrayfun(@(x) {num2str(x)}, 1:num_scales);

    % % Spectral entropy (SpecEn) - Over time
    % case 'SpecEn' 
    %     disp('Computing spectral entropy...')
    %     progressbar('Channels')
    %     parfor iChan = 1:nchan
    %         % fprintf('Channel %d \n', iChan)
    %         [entropy(iChan,:),te] = pentropy( zscore(EEG.data(chanIdx(iChan),:)), EEG.srate);
    % 
    %         fprintf('   %s: %6.3f \n', EEG.chanlocs(iChan).labels, entropy(iChan,:))
    %         progressbar(iChan/nchan)
    %     end

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

% Visualize / Plot
if vis
    if nChan>1
        if strcmpi(measure, 'ExSEnt')
            escape_plot(HD, chanlocs, 'SampEn of durations', scales);
            escape_plot(HA, chanlocs, 'SampEn of amplitudes', scales);
            escape_plot(HDA, chanlocs, 'Joint SampEn of durations & amplitudes', scales);
        else
            % % Standard plot
            escape_plot(entropy, chanlocs, measure, scales);

            % Time-resolved plot
            % escape_plot(info.mse_time, EEG.chanlocs, 'MSE (time-resolved, 2s/1s)', scales, info.time_sec);
        end
    else
        disp("Sorry, you need more than 1 EEG channel for visualization.")
    end
end

% outputs to export in EEG structure
if strcmpi(measure, 'ExSEnt')
    EEG.escape.(measure).data.HD = HD;
    EEG.escape.(measure).data.HA = HA;
    EEG.escape.(measure).data.HDA = HDA;
else
    EEG.escape.(measure).data = entropy;
end
EEG.escape.(measure).electrode_labels = chanLabels;
EEG.escape.(measure).electrode_locations = chanlocs;
if contains(lower(measure), {'mse' 'mfe' 'rcmfe'})
    EEG.escape.(measure).scales = scales;
end

% ADD PARAMETERS USED IN STRUCTURE OUTPUT 


% Command history (TO FIX)
chanLabels = strjoin(chanlist);
chanLabels = insertBefore(chanLabels," ", "'");
chanLabels = insertAfter(chanLabels," ", "'");
com = sprintf('EEG = escape_compute(''%s'', {''%s''}, %d, %d, %s, %d, %d, %s, %d);', ...
    measure,chanLabels,tau,m,coarsing,num_scales,filt_scales,'[]',vis);

disp('Done computing with Escape! Outputs can be found in the EEG.escape structure.')
fprintf('Time to compute: %.2f minutes. \n', toc(tstart)/60)

