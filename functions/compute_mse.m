function [mse, scales, info] = compute_mse(data, varargin)
% COMPUTE_MSE  Multichannel Multiscale Entropy aligned with:
%   • Azami (2017): SampEn uses r = 0.15 * std after z-scoring (PER CHANNEL, PER SCALE).
%   • Costa (2016/2017): optional 'variance' coarse-graining (volatility MSE).
%   • Kosciessa (2020): narrowband annuli + filter-skip to avoid broadband confounds.
%
%   [mse, scales, info] = compute_mse(data, 'Fs', fs, 'm', 2, 'Tau', 1, 'nScales', 20, ...)
%
% REQUIRED
%   data   [nChan x nSamples] numeric matrix  OR  EEGLAB struct with fields .data, .srate
%
% NAME–VALUE PAIRS (all optional unless marked required)
%   'Fs'                (required unless data is EEGLAB struct) sampling frequency in Hz
%   'm'                 embedding dimension for SampEn (default = 2)
%   'Tau'               time delay for SampEn (default = 1)
%   'nScales'           number of scales (default = 20; truncated to keep >= MinSamplesPerBin)
%   'CoarseType'        'sd' (default, Azami) | 'variance' (Costa 2017) | 'mean' | 'median' (robust)
%   'FilterMode'        'narrowband' (default, annuli) | 'lowpass' | 'highpass' | 'none'
%   'FilterDesign'      'fir' (default, recommended for narrowbands) | 'iir'
%   'TransWidth'        FIR transition width in Hz (default auto: max(0.5, 0.15*passband_width))
%   'FIRMaxOrder'       cap on FIR order (default = 2000)
%   'IIROrder'          Butterworth section order if 'iir' (default = 6; avoid for narrowbands)
%   'PadLen'            reflection padding length (samples). Auto: min(max(100, 9*effOrder), floor((T-1)/2))
%   'MinSamplesPerBin'  minimum coarse bins per scale (default = max(4, m+1))
%   'Parallel'          logical true/false (default = false). When true, parfor over channels is used.
%                        (compute_SampEn is always called with 'Parallel','off' to avoid double-parallel)
%   'Progress'          logical true/false (default = true).
%                        • Parallel==true  & Progress==true → text only
%                        • Parallel==false & Progress==true → text + waitbar (fallback to text if headless)
%
% IMPORTANT
%   • Do NOT pass r here. compute_SampEn must z-score internally and use r = 0.15*std(data_z) PER CHANNEL.
%
% OUTPUTS
%   mse     [nChan x nScales_kept] multiscale sample entropy values
%   scales  1 x nScales_kept cell array with [lo hi] Hz per scale
%   info    struct: options, coarse-bin counts (cg_len), elapsed time, fs
%
% EXAMPLES
%   % A) EEGLAB struct input (Fs auto from EEG.srate), parallel over channels
%   [mse, scales] = compute_mse(EEG, 'nScales', 18, 'Parallel', true, 'Progress', true);
%
%   % B) Numeric matrix input + explicit Fs, variance coarse, lowpass
%   [mse, scales] = compute_mse(EEG.data, 'Fs', EEG.srate, 'CoarseType','variance', 'FilterMode','lowpass');
%
%   % C) Robust median coarse-graining, no filtering
%   [mse, scales] = compute_mse(EEG, 'FilterMode','none', 'CoarseType','median', 'nScales', 25);
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

tStart = tic;

% -------- Parse options (accept Fs later, may come from EEG struct)
p = inputParser;
p.addParameter('Fs', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Tau', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('nScales', 20, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('CoarseType', 'sd');
p.addParameter('FilterMode', 'narrowband');    % Kosciessa default
p.addParameter('FilterDesign', 'fir');         % FIR default for narrowbands
p.addParameter('TransWidth', []);              % Hz; auto if empty
p.addParameter('FIRMaxOrder', 2000);
p.addParameter('IIROrder', 6);
p.addParameter('PadLen', []);                  % auto if empty
p.addParameter('MinSamplesPerBin', []);
p.addParameter('Parallel', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,  @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
o = p.Results;

% -------- Accept EEGLAB struct or numeric
[dat, fs_from_struct] = coerce_data_matrix(data);
if isempty(o.Fs)
    if isempty(fs_from_struct)
        error('Fs must be provided (numeric data) or available as EEG.srate (struct input).');
    else
        o.Fs = fs_from_struct;
    end
end
[nChan, nSamp] = size(dat);
if isempty(o.MinSamplesPerBin), o.MinSamplesPerBin = max(4, o.m+1); end

% -------- CoarseType token (supports robust 'median')
ct = parse_coarse_type(o.CoarseType);

% -------- Scales & bands (simple Kosciessa mapping)
% narrowband annuli: for scale s, passband = [nf/(s+1), nf/s]
nf = o.Fs/2;
makeBand = @(s,mode) switch_lower_mode(mode, nf, s);

% ---- Sanitize scales
S = o.nScales;
if ischar(S) || isstring(S), S = str2double(S); end
if ~isscalar(S) || ~isfinite(S) || S < 1, S = 1; end
S = floor(S);

maxS = floor(nSamp / o.MinSamplesPerBin);
if S > maxS
    warning('Reducing nScales from %d to %d to keep >=%d coarse bins.', S, maxS, o.MinSamplesPerBin);
    S = maxS;
end

% If still zero (e.g., tiny nSamp), fallback to S=1 and disable filtering
if S < 1
    warning(['All requested scales would have < %d coarse bins with current settings. ' ...
             'Falling back to nScales=1 with FilterMode=''none'' (statistical coarse-grain). ' ...
             'Consider: longer epochs, fewer scales, FilterMode=''lowpass'' or lower MinSamplesPerBin.'], ...
             o.MinSamplesPerBin);
    S = 1;
    o.FilterMode = 'none';
end

% Precompute bands
scales = cell(1,S);
for s = 1:S
    bw = makeBand(s, o.FilterMode);
    scales{s} = (bw(2) > bw(1)) * bw + ~(bw(2) > bw(1)) * [NaN NaN]; % keep NaN if invalid
    if o.Progress
        if strcmpi(o.FilterMode,'none')
            fprintf('Scale %2d: FULL [0 %.3f] Hz\n', s, nf);
        else
            fprintf('Scale %2d: PASS [%.3f %.3f] Hz\n', s, scales{s}(1), scales{s}(2));
        end
    end
end

% -------- Pre-allocate & progress
mse    = nan(nChan, S);
cg_len = zeros(1, S);

useWB = (o.Progress && ~o.Parallel && usejava('desktop'));
hWB = [];
if useWB
    try, hWB = waitbar(0,'Computing MSE...','Name','compute_mse'); catch, hWB=[]; end
end
if o.Progress
    fprintf('MSE: %d scales | Filter=%s/%s | Coarse=%s | Parallel(ch)=%d\n', ...
        S, o.FilterMode, o.FilterDesign, ct, o.Parallel);
end
tick = @(k) sprintf('Scales %d/%d (%.0f%%)', k, S, 100*k/max(1,S));

% -------- Scales loop (serial), channels (parfor if requested)
for s=1:S
    % 1) Filtering (FIR default; zero-phase with reflection padding)
    Y = dat;
    band = scales{s};
    if ~strcmpi(o.FilterMode,'none') && ~any(isnan(band))
        switch lower(o.FilterDesign)
            case 'fir'
                Y = zp_fir(Y, o.Fs, band(1), band(2), nf, o.TransWidth, o.FIRMaxOrder, o.PadLen);
            otherwise
                Y = zp_iir(Y, band(1), band(2), nf, o.IIROrder, o.PadLen);
        end
    end

    % 2) Coarse-grain: filter-skip when filtering; else statistical coarse-grain
    if ~strcmpi(o.FilterMode,'none') && ~any(isnan(band))
        idx = 1:s:floor(size(Y,2)/s)*s;
        CG  = Y(:, idx);
        nBins = size(CG,2);
    else
        [CG, nBins] = coarsegrain_stat(Y, s, ct);
    end
    cg_len(s) = nBins;

    % 3) SampEn per channel (SampEn: z-score internally, r=0.15*std, Parallel='off')
    if nBins >= max(4, o.m+1)
        mse_s = nan(nChan,1);
        if o.Parallel && ~isempty(ver('parallel'))
            if o.Progress, fprintf('  [par] computing SampEn over %d channel(s) @ scale %d\n', nChan, s); end
            parfor ch=1:nChan
                cg = CG(ch,:);
                mse_s(ch) = compute_SampEn(cg, 'm', o.m, 'tau', o.Tau, 'Parallel', false, 'Progress', false);
            end
        else
            if o.Progress, fprintf('  [ser] computing SampEn over %d channel(s) @ scale %d\n', nChan, s); end
            for ch=1:nChan
                cg = CG(ch,:);
                mse_s(ch) = compute_SampEn(cg, 'm', o.m, 'tau', o.Tau, 'Parallel', false, 'Progress', false);
            end
        end
        mse(:,s) = mse_s;
    else
        if o.Progress
            fprintf('  [skip] scale %d: nBins=%d (< Min=%d)\n', s, nBins, o.MinSamplesPerBin);
        end
    end

    if useWB && ~isempty(hWB) && isvalid(hWB)
        try, waitbar(s/S, hWB, tick(s)); end
    elseif o.Progress
        fprintf(' - scale %d/%d\n', s, S);
    end
end

% -------- Drop NaN-only scales (but never return 0 columns)
nanScale = all(isnan(mse),1);
if any(nanScale) && sum(~nanScale) >= 1
    mse(:,nanScale) = [];
    scales(nanScale) = [];
    cg_len(nanScale) = [];
elseif all(nanScale) && ~isempty(nanScale)
    keep = find(nanScale, 1, 'first');
    drop = setdiff(1:numel(nanScale), keep);
    mse(:,drop) = [];
    scales(drop) = [];
    cg_len(drop) = [];
end

% -------- Close waitbar & info
if ~isempty(hWB) && isvalid(hWB)
    try, close(hWB); end
end
info = struct('FilterMode',o.FilterMode,'FilterDesign',o.FilterDesign,'TransWidth',o.TransWidth, ...
              'FIRMaxOrder',o.FIRMaxOrder,'IIROrder',o.IIROrder,'PadLen',o.PadLen, ...
              'MinSamplesPerBin',o.MinSamplesPerBin,'Parallel',o.Parallel, ...
              'fs',o.Fs,'cg_len',cg_len,'elapsed',toc(tStart));

end % === compute_mse ===


% ======================= LOCAL HELPERS =======================

function [X, fs] = coerce_data_matrix(data)
% Accept numeric [nChan x nSamples] or EEGLAB EEG struct.
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
% Parse coarse-type; support numeric aliases and robust 'median'.
ct = 'sd';
if isempty(token), return; end
if isnumeric(token)
    switch round(token)
        case 0, ct='mean';
        case 1, ct='sd';
        case 2, ct='variance';
        otherwise, ct='sd';
    end
else
    t = lower(regexprep(strtrim(char(token)),'[^a-z]',''));
    if any(strcmp(t,{'sd','sigma','std','standarddeviation'})), ct='sd';
    elseif any(strcmp(t,{'variance','var','sigma2'})),         ct='variance';
    elseif any(strcmp(t,{'mean','avg','average','mu'})),       ct='mean';
    elseif any(strcmp(t,{'median','med'})),                    ct='median';
    else, ct='sd';
    end
end
end

function band = switch_lower_mode(mode, nf, s)
% Map scale s -> [lo hi] Hz band per FilterMode (Kosciessa-simple).
mode = lower(mode);
switch mode
    case 'narrowband'  % annuli
        band = [max(0, nf/(s+1)), min(nf, nf/s)];
    case 'lowpass'
        band = [0, min(nf, nf/s)];
    case 'highpass'
        band = [max(0, nf/(s+1)), nf];
    otherwise % 'none'
        band = [0, nf];
end
end

function Y = zp_fir(Y, fs, lo, hi, nf, transHz, maxOrder, padLen)
% Zero-phase linear-phase FIR (Hamming), reflection padded.
% Auto order from transition width; capped by maxOrder for speed.
if hi <= lo || lo < 0 || hi > nf, return; end
pbw = max(hi-lo, eps);
if isempty(transHz)
    transHz = max(0.5, 0.15*pbw); % >=0.5 Hz or 15% of passband width
end
tw = max(transHz, fs/size(Y,2)); % guard: at least ~1 bin in transition
ord = min(maxOrder, max(10, ceil(3.3*fs/tw)));  % Hamming ~60 dB

% Build FIR
if lo<=0 && hi>0
    Wn = hi/(fs/2); b  = fir1(ord, Wn, 'low',  hamming(ord+1), 'scale');
elseif hi>=nf && lo>0
    Wn = lo/(fs/2); b  = fir1(ord, Wn, 'high', hamming(ord+1), 'scale');
else
    Wn = [lo hi]/(fs/2); b = fir1(ord, Wn, 'bandpass', hamming(ord+1), 'scale');
end

% Reflection padding + zero-phase filtering
T = size(Y,2);
if isempty(padLen)
    pad = min(max(100, 9*ord), floor((T-1)/2));
else
    pad = min(padLen, floor((T-1)/2));
end
if pad>0
    Ypad = [fliplr(Y(:,1:pad)) , Y , fliplr(Y(:,end-pad+1:end))];
    Ypad = filtfilt(b, 1, Ypad.').';
    Y = Ypad(:, pad+1:end-pad);
else
    Y = filtfilt(b, 1, Y.').';
end
end

function Y = zp_iir(Y, lo, hi, nf, ord, padLen)
% Zero-phase IIR Butterworth (SOS if available) with reflection padding
if hi <= lo || lo < 0 || lo > nf, return; end %#ok<*ISMT>
Wn = [lo hi]/nf; ftype='bandpass';
if lo<=0 && hi>0,  Wn = hi/nf;  ftype='low'; end
if hi>=nf && lo>0, Wn = lo/nf;  ftype='high'; end
if any(Wn<=0) || any(Wn>=1), return; end
useSOS = exist('sosfiltfilt','file')==2;
if useSOS, [z,p,k] = butter(ord, Wn, ftype); sos = zp2sos(z,p,k);
else,      [b,a]    = butter(ord, Wn, ftype); end
T = size(Y,2);
if isempty(padLen)
    effOrder = ord*2;
    pad = min(max(100, 9*effOrder), floor((T-1)/2));
else
    pad = min(padLen, floor((T-1)/2));
end
if pad>0
    Ypad = [fliplr(Y(:,1:pad)) , Y , fliplr(Y(:,end-pad+1:end))];
    if useSOS, Ypad = sosfiltfilt(sos, Ypad.').'; else, Ypad = filtfilt(b,a, Ypad.').'; end
    Y = Ypad(:, pad+1:end-pad);
else
    if useSOS, Y = sosfiltfilt(sos, Y.').'; else, Y = filtfilt(b,a, Y.').'; end
end
end

function [CG, nBins] = coarsegrain_stat(Y, s, ct)
% Statistical coarse-graining by factor s; supports 'sd'|'variance'|'mean'|'median'
nChan = size(Y,1);
nFull = floor(size(Y,2)/s)*s;
if nFull == 0, CG = nan(nChan,0); nBins = 0; return; end
nBins = nFull / s;
Z = reshape(Y(:,1:nFull), nChan, s, nBins);
switch ct
    case 'sd'
        CG = squeeze(std(Z, 0, 2, 'omitnan'));
    case 'variance'
        CG = squeeze(var(Z, 0, 2, 'omitnan'));
    case 'mean'
        CG = squeeze(mean(Z, 2, 'omitnan'));
    case 'median'
        CG = squeeze(median(Z, 2, 'omitnan'));
    otherwise
        CG = squeeze(std(Z, 0, 2, 'omitnan')); % safe fallback
end
if isvector(CG), CG = CG(:)'; end
end
