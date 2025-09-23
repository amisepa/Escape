function [mse, scales, info] = compute_mMSE(data, varargin)
% compute_mse  Multichannel Multiscale Entropy (MSE) with FT-style filtskip,
%              time-resolved output, and multiple coarse-graining statistics.
%
% SUMMARY
% -------
% • Two processing modes
%   (A) Statistical coarse-graining (filter_mode='none'): Costa-like block
%       reduction using the chosen statistic ('mean' | 'std' | 'variance' | 'median').
%   (B) Spectral coarse-graining (filter_mode='narrowband'|'lowpass'|'highpass'):
%       FieldTrip/Kosciessa-style **filtskip** — filter to a scale-specific band,
%       decimate by scale s, and **average SampEn across all s start phases**.
%
% Narrowband behavior (FieldTrip-like):
%   • Scale 1: high-pass only at ~Nyq/2 (widened downward), NO skipping.
%   • Scales >1: bandpass on annuli [Nyq/(s+1), Nyq/s] widened by ±NBWiden,
%     then skip every s samples and average across all s start phases.
%   • Auto-switch to IIR for narrow bands; else zero-phase FIR with robust padding.
%
% Scales are auto-dropped when infeasible (too few coarse bins, or FIR band unrealizable).
%
%
% ADVANTAGES vs. traditional Costa method (2002)
% ---------------------------
% 1) Reduced large-scale bias: Costa shortens the series to N/s. Here we enforce
%    sensible bin counts and (for filtskip) average across s phase starts,
%    yielding smoother, more stable curves at large scales.
% 2) Richer coarse statistics: beyond 'mean', support 'std'/'variance'/'median'
%    (volatility-MSE variants useful for amplitude/variance dynamics).
% 3) Time-resolved analysis: sliding-window MSE at each scale (also in filtskip).
% 4) Narrowband MSE (Kosciessa et al.): isolate band-limited entropy via scale
%    annuli with ±NBWiden; scale 1 uses HP-only.
%
%
% DIFFERENCES relative to CMSE/RCMSE (Azami)
% ------------------------------------------
% • Filtskip mode here is **composite across phases** by averaging *entropy values*
%   from all s start phases (CMSE-like). **RCMSE** pools the *match counts*
%   (ϕ_m, ϕ_{m+1}) across phases before taking −ln, further reducing variance and
%   avoiding undefined values. This function does **not** implement that refined pooling.
% • Tolerance r: this function forwards 'r' to SampEn as a **fraction of the SD of
%   each (coarse-grained or filtered) series**. This often flattens white-noise MSE.
%   To reproduce classic Costa/RCMSE "white decreases with scale", use a SampEn
%   variant that fixes r to the original (scale-1) SD.
%
%
% INPUTS (name–value)
% -------------------
%   'Fs'              : sampling rate in Hz (required unless EEG.srate present)
%   'm'               : embedding dimension (default 2, ≥1)
%   'tau'             : lag (default 1, ≥1)
%   'r'               : tolerance **fraction of SD** for SampEn (default 0.15)
%   'num_scales'      : requested scales; infeasible scales are dropped (default 20)
%   'coarsing'        : 'mean' | 'std' | 'variance' | 'median' (default 'std')
%   'filter_mode'     : 'none' | 'narrowband' | 'lowpass' | 'highpass' (default 'none')
%   'FilterDesign'    : 'fir' | 'iir' (default 'fir'; auto IIR for very narrow bands)
%   'TransWidth'      : FIR transition width in Hz (default [], auto heuristic)
%   'FIRMaxOrder'     : max FIR order (default 2000)
%   'IIROrder'        : IIR (Butterworth) order (default 10)
%   'PadLen'          : reflection padding for filtfilt (samples; default [], auto)
%   'MinSamplesPerBin': minimum coarse bins per scale (default max(4, m+1))
%   'Parallel'        : parallelize across channels (default true)
%   'Progress'        : console/waitbar feedback (default true)
%   'NBWiden'         : ±fraction widening for filtskip bands (default 0.05; 0 ≤ < 0.25)
%   'TimeWin'         : time-window length in seconds for time-resolved MSE (default [])
%   'TOI'             : vector of window centers (s); if empty, auto grid (default [])
%   'TimeStep'        : step between centers (s); default = TimeWin/2 (when TimeWin set)
%
% Accepted data:
%   data : double [nChan×nSamples] or EEGLAB EEG struct (.data, .srate)
%
%
% OUTPUTS
% -------
%   mse      : double [nChan×nScales_kept] — MSE (SampEn) per channel & scale
%   scales   : 1×nScales_kept cellstr — numeric scales (statistical) or band labels (filtskip)
%   info     : struct with fields:
%              • filter_mode, FilterDesign, TransWidth, FIRMaxOrder, IIROrder, PadLen
%              • MinSamplesPerBin, Parallel, fs, elapsed
%              • cg_len (1×nScales_kept): effective coarse length per scale
%              • NBWiden, TimeWin, TOI, TimeStep
%              • mse_time (optional): double [nChan×nScales_kept×nTime] — time-resolved MSE
%              • time_sec (optional): 1×nTime — window center times (s)
%
%
% USAGE EXAMPLES
% --------------
% 1) Costa-like MSE (statistical path), mean coarse-grain
%    [mse, scales, info] = compute_mse(X, 'Fs', 256, 'm', 2, 'tau', 1, ...
%        'r', 0.15, 'num_scales', 30, 'coarsing', 'mean', ...
%        'filter_mode','none', 'Parallel', true, 'Progress', true);
%
% 2) Volatility-MSE (std), statistical path
%    [mse, scales] = compute_mse(X, 'Fs', 512, 'coarsing','std', ...
%        'num_scales', 25, 'filter_mode','none');
%
% 3) Narrowband filtskip (Kosciessa-style), FIR when feasible
%    [mse, scales] = compute_mse(X, 'Fs', 500, 'filter_mode','narrowband', ...
%        'NBWiden', 0.05, 'FilterDesign','fir', 'num_scales', 20, 'coarsing','mean');
%
% 4) Time-resolved narrowband MSE, 2-s windows every 0.5-s
%    [mse, scales, info] = compute_mse(X, 'Fs', 256, 'filter_mode','narrowband', ...
%        'TimeWin', 2.0, 'TimeStep', 0.5, 'num_scales', 15);
%    % → info.mse_time: [nChan×nScales×nTime], info.time_sec: centers (s)
%
% 5) EEGLAB struct input with auto Fs
%    [mse, scales] = compute_mse(EEG, 'num_scales', 20, 'filter_mode','none');
%
%
% Notes
% -----
% • 'std'/'variance' at scale 1 can be numerically trivial; guards keep unstable
%   cases from propagating. Many analyses start volatility-MSE at scale ≥ 2.
% • In filtskip mode, scale 1 is HP-only; scales >1 are annuli widened by ±NBWiden.
%   Extremely narrow bands auto-switch to IIR; unrealizable FIR bands are dropped.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------


tStart = tic;

% -------- Parse options
p = inputParser;
p.addParameter('Fs', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('r', .15, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('num_scales', 20, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('coarsing', 'std');
p.addParameter('filter_mode', 'none');      % 'narrowband' | 'lowpass' | 'highpass' | 'none'
p.addParameter('FilterDesign', 'fir');            % 'fir' | 'iir'
p.addParameter('TransWidth', []);                 % FIR desired transition (Hz); auto if empty
p.addParameter('FIRMaxOrder', 2000);
p.addParameter('IIROrder', 10);                   % steeper for BP (closer to FT)
p.addParameter('PadLen', []);                     % reflection padding length; auto if empty
p.addParameter('MinSamplesPerBin', []);
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,  @(x) islogical(x) && isscalar(x));
p.addParameter('NBWiden', 0.05, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<0.25); % ±5% widen
p.addParameter('TimeWin', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('TOI', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('TimeStep', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.parse(varargin{:});
o = p.Results;

% -------- Accept EEGLAB struct or numeric
[dat, fs_from_struct] = coerce_data_matrix(data);
if isempty(o.Fs)
    if isempty(fs_from_struct)
        error('Fs must be provided (numeric) or present in EEG.srate (struct input).');
    else
        o.Fs = fs_from_struct;
    end
end
[nChan, nSamp] = size(dat);
if isempty(o.MinSamplesPerBin), o.MinSamplesPerBin = max(4, o.m+1); end
ct = parse_coarse_type(o.coarsing);

% -------- Time grid (if requested)
Tsec   = nSamp / o.Fs;
doTime = ~isempty(o.TimeWin);
if doTime
    if isempty(o.TOI)
        step = iff(isempty(o.TimeStep), o.TimeWin/2, o.TimeStep);
        t0 = o.TimeWin/2; t1 = Tsec - o.TimeWin/2;
        if t1 < t0
            warning('TimeWin is larger than the data duration; disabling time-resolved output.');
            doTime = false;
        else
            tCenters = t0:step:t1;
        end
    else
        tCenters = o.TOI(:).';
        half = o.TimeWin/2;
        tCenters = tCenters(tCenters>=half & tCenters<=Tsec-half);
        if isempty(tCenters), doTime = false; end
    end
else
    tCenters = [];
end
nTOI = numel(tCenters);

% -------- Determine usable number of scales
S    = max(1, floor(o.num_scales));
maxS = floor(nSamp / o.MinSamplesPerBin);
if S > maxS
    warning('Reducing num_scales from %d to %d to keep >=%d coarse bins.', S, maxS, o.MinSamplesPerBin);
    S = maxS;
end
if S < 1
    S = 1; o.filter_mode = 'none';
end

% -------- FT-style bands (±NBWiden); s=1 HP-only
nyq   = o.Fs/2;
bands = cell(1,S);
scales= cell(1,S);
if strcmpi(o.filter_mode,'none')
    for s=1:S, bands{s} = [NaN NaN]; scales{s} = num2str(s); end
else
    for s=1:S
        lo_raw = nyq/(s+1);
        hi_raw = nyq/s;
        if s==1
            lo_eff   = max(0, lo_raw*(1 - o.NBWiden)); % widen downward
            bands{s} = [lo_eff, Inf];                  % HP-only flag
            scales{s}= sprintf('HP>%.3f Hz', lo_eff);
        else
            lo_eff   = max(0, lo_raw*(1 - o.NBWiden));
            hi_eff   = min(nyq*(1-1e-6), hi_raw*(1 + o.NBWiden));
            if hi_eff <= lo_eff
                bands{s} = [NaN NaN]; scales{s}='invalid';
            else
                bands{s} = [lo_eff, hi_eff];
                scales{s}= sprintf('[%.3f %.3f] Hz', lo_eff, hi_eff);
            end
        end
        if o.Progress && all(isfinite(bands{s}))
            fprintf('Scale %2d band: %s\n', s, scales{s});
        end
    end
end

% -------- Pre-allocate & progress
mse      = nan(nChan, S);
cg_len   = zeros(1, S);
mse_time = []; if doTime, mse_time = nan(nChan, S, nTOI); end

useWB = (o.Progress && ~o.Parallel && usejava('desktop'));
hWB = [];
if useWB
    try
        hWB = waitbar(0,'Computing MSE...','Name','compute_mse progress bar');
        waitbar(0, hWB, sprintf('Scale 1/%d • Channel 0/%d', S, nChan));
    catch, hWB = []; end
end
if o.Progress
    fprintf('MSE: %d scales | Filter=%s/%s | Coarse=%s | Parallel(ch)=%d | TimeWin=%s | nTOI=%d\n', ...
        S, o.filter_mode, o.FilterDesign, ct, o.Parallel, tern(doTime, sprintf('%.3gs', o.TimeWin), 'off'), nTOI);
end

% -------- Fill missing values pre-filter
for ch = 1:nChan
    x = dat(ch,:); if any(~isfinite(x)), dat(ch,:) = fillmissing(x,'linear','EndValues','nearest'); end
end

% -------- Scales loop
for s=1:S
    feff = o.Fs / s;
    Y    = dat;
    band = bands{s};

    % --- Statistical path (no filtering)
    if strcmpi(o.filter_mode,'none') || any(isnan(band))
        [CG, nBins] = coarsegrain_stat(Y, s, ct);
        cg_len(s)   = nBins;
        minNeeded   = max([o.MinSamplesPerBin, o.m+1, 100]);  % Grandy stability floor
        if nBins < minNeeded
            if o.Progress, fprintf('  [drop] scale %d: nBins=%d (<%d)\n', s, nBins, minNeeded); end
            continue
        end

        vals = nan(nChan,1);
        if o.Parallel && ~isempty(ver('parallel'))
            parfor ch = 1:nChan
                vals(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                                          'r', o.r, 'Parallel', false, 'Progress', false);
            end
        else
            for ch = 1:nChan
                vals(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                                          'r', o.r, 'Parallel', false, 'Progress', false);
            end
        end
        mse(:,s) = vals;

        if doTime
            mse_time = time_windows(CG, s, feff, o, mse_time, nTOI, tCenters, hWB, S, nChan);
        end
        continue
    end

    % --- Filtering (FT-like); auto IIR for very narrow bands
    isHPonly = isinf(band(2));
    if strcmpi(o.FilterDesign,'iir')
        Y = ft_iir(Y, band, o.Fs, o.IIROrder, o.PadLen);
    else
        % Heuristic: prefer IIR when band is extremely narrow in Hz
        if isHPonly
            pbw = nyq - band(1);
        else
            pbw = band(2) - band(1);
        end
        veryNarrow = (pbw / o.Fs) < 0.01;   % <1% of Fs
        if veryNarrow
            Y = ft_iir(Y, band, o.Fs, o.IIROrder, o.PadLen);
        else
            Y = fir_zero_phase(Y, band, o.Fs, o.TransWidth, o.FIRMaxOrder, o.PadLen);
        end
    end

    % --- Filtskip/coarse-grain and feasibility drop
    if isHPonly
        nStarts   = 1;
        nBins_eff = size(Y,2);
    else
        nStarts   = s;
        nBins_eff = floor(size(Y,2)/s);
    end
    cg_len(s) = nBins_eff;

    % drop guard: too few coarse bins OR FIR unrealizable (if FIR path)
    % tw_guess = feval(@()0); % dummy init
    if strcmpi(o.FilterDesign,'fir')
        if isHPonly
            pbw_guard = nyq - band(1);
        else
            pbw_guard = band(2) - band(1);
        end
        if isempty(o.TransWidth)
            tw_guess = max(0.15*max(pbw_guard,eps), o.Fs/size(Y,2));
        else
            tw_guess = max(o.TransWidth, o.Fs/size(Y,2));
        end
        if (~isHPonly) && (pbw_guard < 2*tw_guess)
            if o.Progress, fprintf('  [drop] scale %d: FIR band too narrow (pbw=%.4gHz < 2*tw≈%.4gHz)\n', s, pbw_guard, tw_guess); end
            continue
        end
    end
    if nBins_eff < max(4, o.m+1)
        if o.Progress, fprintf('  [drop] scale %d: nBins=%d (<%d)\n', s, nBins_eff, max(4,o.m+1)); end
        continue
    end

    % --- Whole-epoch SampEn across starts
    vals = nan(nChan, nStarts);
    if o.Parallel && ~isempty(ver('parallel'))
        parfor is = 1:nStarts
            if nStarts==1
                CG = Y;
            else
                last = is+(nBins_eff-1)*s;
                idx  = is:s:last;
                CG   = Y(:, idx);
            end
            % vals(:,is) = sampen_block(CG, o.m, o.tau, o.r, false);
            tmp = nan(nChan,1);
            for ch = 1:nChan
                tmp(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                                         'r', o.r, 'Parallel', false, 'Progress', false);
            end
            vals(:,is) = tmp;
        end
    else
        for is = 1:nStarts
            if nStarts==1
                CG = Y;
            else
                last = is+(nBins_eff-1)*s;
                idx  = is:s:last;
                CG   = Y(:, idx);
            end

            tmp = nan(nChan,1);
            for ch = 1:nChan
                tmp(ch) = compute_SampEn(CG(ch,:), 'm', o.m, 'tau', o.tau, ...
                                         'r', o.r, 'Parallel', false, 'Progress', false);
            end
            vals(:,is) = tmp;

            if ~doTime && ~isempty(hWB) && isvalid(hWB)
                try waitbar(((s-1)*nChan + nChan)/(S*nChan), hWB); catch; end
            end
        end
    end
    mse(:,s) = mean(vals,2,'omitnan');

    % --- Time-resolved (average across starts)
    if doTime
        if nStarts==1
            feff_eff = o.Fs;     % HP-only, no skipping
        else
            feff_eff = feff;
        end
        mse_time = time_windows_filtskip(Y, s, nStarts, nBins_eff, feff_eff, o, mse_time, nTOI, tCenters, hWB, S, nChan);
    end
end

% Drop NaN-only scales (but never return 0 columns)
nanScale = all(isnan(mse),1);
if any(nanScale) && sum(~nanScale) >= 1
    mse(:,nanScale) = []; scales(nanScale) = []; cg_len(nanScale) = [];
    if doTime, mse_time(:,nanScale,:) = []; end
elseif all(nanScale) && ~isempty(nanScale)
    keep = find(nanScale, 1, 'first');
    drop = setdiff(1:numel(nanScale), keep);
    mse(:,drop) = []; scales(drop) = []; cg_len(drop) = [];
    if doTime, mse_time(:,drop,:) = []; end
end

% Close waitbar & info
if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch; end, end
info = struct('filter_mode',o.filter_mode,'FilterDesign',o.FilterDesign,'TransWidth',o.TransWidth, ...
    'FIRMaxOrder',o.FIRMaxOrder,'IIROrder',o.IIROrder,'PadLen',o.PadLen, ...
    'MinSamplesPerBin',o.MinSamplesPerBin,'Parallel',o.Parallel, ...
    'fs',o.Fs,'cg_len',cg_len,'elapsed',toc(tStart), ...
    'NBWiden',o.NBWiden,'TimeWin',o.TimeWin,'TOI',o.TOI,'TimeStep',o.TimeStep);

if doTime
    info.mse_time = mse_time;          % chan x scale x time
    if ~isfield(info,'time_sec') || isempty(info.time_sec)
        info.time_sec = tCenters;      % if no trimming occurred
    end
end

end % compute_mse


%% ======================= LOCAL HELPERS =======================

function y = iff(cond, a, b)
if cond, y=a; else, y=b; end
end
function y = tern(cond, a, b)
if cond, y=a; else, y=b; end
end

function [X, fs] = coerce_data_matrix(data)
fs = [];
if isnumeric(data)
    X = double(data);
elseif isstruct(data) && isfield(data,'data')
    X = double(data.data);
    if isfield(data,'srate') && ~isempty(data.srate), fs = double(data.srate); end
else
    error('First argument must be numeric [nChan x nSamples] or EEGLAB EEG struct with field .data');
end
if ~isreal(X)
    warning('Data is complex; using real part.'); X = real(X);
end
end

function ct = parse_coarse_type(token)
ct = 'std';
if isempty(token), return; end
if isnumeric(token)
    switch round(token)
        case 0, ct='mean';
        case 1, ct='std';
        case 2, ct='variance';
        otherwise, ct='std';
    end
else
    t = lower(regexprep(strtrim(char(token)),'[^a-z]',''));
    if any(strcmp(t,{'std','sigma','std','standarddeviation'})), ct='std';
    elseif any(strcmp(t,{'variance','var','sigma2'})),         ct='variance';
    elseif any(strcmp(t,{'mean','avg','average','mu'})),       ct='mean';
    elseif any(strcmp(t,{'median','med'})),                    ct='median';
    else, ct='std';
    end
end
end

function [CG, nBins] = coarsegrain_stat(Y, s, ct)
nChan = size(Y,1);
nFull = floor(size(Y,2)/s)*s;
if nFull == 0, CG = nan(nChan,0); nBins = 0; return; end
nBins = nFull / s;
Z = reshape(Y(:,1:nFull), nChan, s, nBins);
switch ct
    case 'std',       CG = squeeze(std(Z, 0, 2, 'omitnan'));
    case 'variance', CG = squeeze(var(Z, 0, 2, 'omitnan'));
    case 'mean',     CG = squeeze(mean(Z, 2, 'omitnan'));
    case 'median',   CG = squeeze(median(Z, 2, 'omitnan'));
    otherwise,       CG = squeeze(std(Z, 0, 2, 'omitnan'));
end
if isvector(CG), CG = CG(:)'; end
end


function mse_time = time_windows(CG, s, feff, o, mse_time, nTOI, tCenters, hWB, S, nChan)
W = max(1, round(o.TimeWin * feff));
minW = max([o.m+1, 100]); % 'none' path uses 100-pt floor
if W < minW
    if o.Progress, fprintf('  [skip] scale %d: TimeWin too short (W=%d < %d)\n', s, W, minW); end
    return
end
centers_bins = round(tCenters * feff) + 1;
nBins = size(CG,2);
starts = centers_bins - floor(W/2);
stops  = starts + W - 1;
keep   = starts>=1 & stops<=nBins;
if any(keep)
    baseIdx = find(keep,1,'first');
    for ti = 1:sum(keep)
        seg = CG(:, starts(baseIdx+ti-1):stops(baseIdx+ti-1));
        for ch=1:nChan
            mse_time(ch,s, baseIdx+ti-1) = compute_SampEn(seg(ch,:), 'm', o.m, ...
                'tau', o.tau, 'r', o.r, 'Parallel', false, 'Progress', false);
            if ~isempty(hWB) && isvalid(hWB)
                try, waitbar(((s-1)*nChan + ch)/(S*nChan), hWB); end
            end
        end
    end
end
end

function mse_time = time_windows_filtskip(Y, s, nStarts, nBins_eff, feff, o, mse_time, nTOI, tCenters, hWB, S, nChan)
W = max(1, round(o.TimeWin * feff));
if W < max(4, o.m+1)
    if o.Progress, fprintf('  [skip] scale %d: TimeWin too short (W=%d)\n', s, W); end
    return
end
centers_bins = round(tCenters * feff) + 1;
if nStarts==1
    nBinsAvail = size(Y,2);
else
    nBinsAvail = nBins_eff;
end
starts = centers_bins - floor(W/2);
stops  = starts + W - 1;
keep   = starts>=1 & stops<=nBinsAvail;
if any(keep)
    baseIdx = find(keep,1,'first');
    for ti = 1:sum(keep)
        tiAbs = baseIdx + ti - 1;
        acc = zeros(nChan,1);
        for is = 1:nStarts
            if nStarts==1
                idxv = starts(tiAbs):stops(tiAbs);
            else
                idx0 = is + (starts(tiAbs)-1)*s;
                idxv = idx0 : s : idx0 + (W-1)*s;
                idxv(idxv>size(Y,2))=[];
            end
            if numel(idxv) >= max(4, o.m+1)
                seg = Y(:, idxv);
                for ch=1:nChan
                    acc(ch) = acc(ch) + compute_SampEn(seg(ch,:), 'm', o.m, ...
                        'tau', o.tau,'r', o.r, 'Parallel', false, 'Progress', false);
                end
            else
                acc(:) = NaN; break
            end
        end
        mse_time(:,s, tiAbs) = acc ./ nStarts;
        if ~isempty(hWB) && isvalid(hWB)
            try, waitbar(((s-1)*nChan + nChan)/(S*nChan), hWB); end
        end
    end
end
end

%% ---------- Filters (FIR/IIR) ----------

function Y = fir_zero_phase(Y, band, fs, transHz, maxOrder, padLen)
% Robust fir1 design (no explicit window), handles Nyquist parity internally.
nyq = fs/2;
lo = band(1); hi = band(2);  % hi==Inf → HP-only
if isinf(hi)
    Wn = max(eps, lo/nyq);  ftype = 'high';
elseif lo<=0
    Wn = min(0.999999, hi/nyq); ftype = 'low';
else
    Wn = [max(eps, lo/nyq), min(0.999999, hi/nyq)]; ftype='bandpass';
end

% Order from desired transition; let fir1 adjust if needed.
T  = size(Y,2);
if isempty(transHz), transHz = max(0.5, 0.15*max((isinf(hi))* (nyq-lo) + (~isinf(hi))*(hi-lo), eps)); end
tw = max(transHz, fs/T);
ord = min(maxOrder, max(10, ceil(3.3*fs/tw)));
if mod(ord,2)==1, ord = ord + 1; end
if strcmp(ftype,'bandpass') && mod(ord,2)==1, ord = ord+1; end  % avoid Nyquist warning

b = fir1(ord, Wn, ftype, 'scale');   % no explicit window → no length mismatch if ord changes
L = length(b)-1;

% Reflection padding proportional to actual filter length
if isempty(padLen), pad = min(max(100, 9*L), floor((T-1)/2)); else, pad = min(padLen, floor((T-1)/2)); end
if pad>0
    Ypad = [fliplr(Y(:,1:pad)) , Y , fliplr(Y(:,end-pad+1:end))];
    Ypad = filtfilt(b, 1, Ypad.').'; Y = Ypad(:, pad+1:end-pad);
else
    Y = filtfilt(b, 1, Y.').';
end
end

function Y = ft_iir(Y, band, fs, ord, padLen)
nyq = fs/2; lo = band(1); hi = band(2);
if isinf(hi)
    [b,a] = butter(max(3,ceil(ord/2)), max(eps, lo/nyq), 'high');
    Y = iir_filtpad(Y, b, a, padLen);
else
    [bl,al] = butter(ord, min(0.999999, hi/nyq), 'low');
    [bh,ah] = butter(ord, max(eps,    lo/nyq), 'high');
    Y = iir_filtpad(Y, bl, al, padLen);
    Y = iir_filtpad(Y, bh, ah, padLen);
end
end

function Y = iir_filtpad(Y, b, a, padLen)
T = size(Y,2);
effOrder = max(length(a),length(b));
if isempty(padLen), pad = min(max(100, 9*effOrder), floor((T-1)/2)); else, pad = min(padLen, floor((T-1)/2)); end
if pad>0
    Ypad = [fliplr(Y(:,1:pad)) , Y , fliplr(Y(:,end-pad+1:end))];
    Ypad = filtfilt(b,a, Ypad.').'; Y = Ypad(:, pad+1:end-pad);
else
    Y = filtfilt(b,a, Y.').';
end
end
