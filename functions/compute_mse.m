function [mse, scales, info] = compute_mse(X, fs, m, r, tau, coarseType, nScales, filtData, varargin)
% COMPUTE_MSE  Multichannel Multiscale Entropy with scale-wise similarity bounds and filter-defined scales.
%
% Usage:
%   [mse, scales, info] = compute_mse(X, fs, m, r, tau, coarseType, nScales, filtData, ...)
%
% Inputs:
%   X          - Data matrix [nChan x nSamples]. Each row is one channel.
%   fs         - Sampling frequency (Hz).
%   m          - Embedding dimension for sample entropy (commonly m=2).
%   r          - Similarity criterion (as a RATIO). At each scale, epsilon = r * std(coarse-grained series).
%   tau        - Time delay for embedding (usually 1).
%   coarseType - Coarse-graining statistic (case/alias tolerant):
%                   'Mean' / 'mean' / 0
%                   'SD' | 'Standard deviation' | 'standard deviation' | 'sigma' | 'std' / 1   (default, recommended)
%                   'Variance' | 'variance' | 'var' / 2
%   nScales    - Maximum scale factor to compute.
%   filtData   - If true, applies zero-phase filtering per scale (LP/HP/BP/BS) to limit spectral bias.
%
% Name-value pairs:
%   'FilterMode'       : 'lowpass' (default) | 'highpass' | 'bandpass' | 'bandstop' | 'none'
%   'Band'             : [lo hi] Hz for bandpass/bandstop (clamped to [0, fs/2]).
%   'RescaleBounds'    : true (default). Use per-scale epsilon = r * std(cg).
%                        If false, uses global epsilon = r * sd0(ch) (not recommended).
%   'MinSamplesPerBin' : Minimum coarse-grained LENGTH (bins) to allow per scale (default = max(4, m+1)).
%   'FiltOrder'        : IIR order per section (default = 6). SOS zero-phase used when available.
%   'Parallel'         : 'none' (default) | 'scale' | 'channel'. (After vectorization, 'scale' is usually fastest.)
%   'Verbose'          : true/false (prints band + stats per scale and QC skip reasons).
%   'Progress'         : true/false (default = true). Show progress (text or waitbar).
%   'ProgressStyle'    : 'text' (default) | 'waitbar'. Waitbar needs desktop graphics.
%   'BandMode'         : 'fixed' (default) | 'annuli'  (annuli = non-overlapping rings for bandpass/bandstop)
%   'CoarseMethod'     : 'stat' (Mean/SD/Variance) | 'skip' (filt-skip after filtering)
%
% Outputs:
%   mse     - [nChan x nScales_kept] matrix of multiscale sample entropy values (NaN-only scales removed).
%   scales  - 1 x nScales_kept cell array of [lo hi] Hz per scale (after NaN-only removal).
%   info    - Struct with runtime options and QC:
%               .FilterMode, .Band, .RescaleBounds, .MinSamplesPerBin, .FiltOrder, .Parallel, .Verbose, .fs
%               .cg_len  -> 1 x nScales original coarse-grained lengths (bins) per scale
%               .elapsed -> elapsed time (seconds)
%
% Integrity notes:
%   • SampEn uses absolute epsilon **per scale** (epsilon = r * std(coarse series)), preventing variance-driven bias.
%   • Zero-phase SOS Butterworth with reflection padding minimizes phase and edge artifacts.
%   • Filtering is vectorized across channels per scale; helpers are subfunctions (parfor-safe).
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape plugin
% Author: Cedric Cannard
% GNU GPL v2 or later.
% -------------------------------------------------------------------------

% -------- Options --------
p = inputParser;
p.addParameter('FilterMode','lowpass');                % 'lowpass'|'highpass'|'bandpass'|'bandstop'|'none'
p.addParameter('Band',[8 12]);                         % used for bandpass/bandstop
p.addParameter('RescaleBounds', true);                 % epsilon per scale = r * std(cg)
p.addParameter('MinSamplesPerBin', max(4,m+1));        % minimum # coarse bins per scale
p.addParameter('FiltOrder', 6);                        % Butterworth order per section
p.addParameter('Parallel', 'none');                    % 'none'|'scale'|'channel'
p.addParameter('Verbose', false);
p.addParameter('Progress', true);
p.addParameter('ProgressStyle', 'text');               % 'text'|'waitbar'
p.addParameter('BandMode','fixed');                    % 'fixed'|'annuli'
p.addParameter('CoarseMethod','stat');                 % 'stat'|'skip'
p.parse(varargin{:});
o = p.Results;

tStart = tic;

[nChan, nSamp] = size(X);
nf = fs/2;

% -------- CoarseType parsing (default SD) --------
ct = 'sd';
if nargin >= 6 && ~isempty(coarseType)
    if ischar(coarseType) || isstring(coarseType)
        t = lower(char(coarseType)); t = strtrim(t); t = regexprep(t,'[^a-z]','');
        if any(strcmp(t,{'mean','mu'})), ct='mean';
        elseif any(strcmp(t,{'sd','sigma','standarddeviation','std'})), ct='sd';
        elseif any(strcmp(t,{'variance','var','sigma2'})), ct='variance';
        else, warning('Unknown coarseType "%s"; defaulting to SD.', t); ct='sd';
        end
    elseif isnumeric(coarseType) && isscalar(coarseType)
        switch round(coarseType), case 0, ct='mean'; case 1, ct='sd'; case 2, ct='variance';
            otherwise, warning('Unknown numeric coarseType=%g; defaulting to SD.', coarseType); ct='sd';
        end
    else
        warning('Unrecognized coarseType; defaulting to SD.'); ct='sd';
    end
end

% -------- Sanitize nScales (must be finite integer ≥1); guard by MinSamplesPerBin --------
if ischar(nScales) || isstring(nScales), nScales = str2double(nScales); end
if ~isscalar(nScales) || ~isfinite(nScales), error('nScales must be a finite scalar.'); end
nScales = floor(double(nScales));
maxScales = floor(nSamp / o.MinSamplesPerBin);
if nScales > maxScales
    warning('Reducing nScales from %d to %d (>= %d bins/scale).', nScales, maxScales, o.MinSamplesPerBin);
    nScales = maxScales;
end
nScales = max(1,nScales);
if strcmpi(o.Parallel,'scale') && nScales <= 1, o.Parallel = 'none'; end

% -------- Warn if skip without filtering --------
if strcmpi(o.CoarseMethod,'skip') && ~filtData
    warning('CoarseMethod="skip" without filtering is not recommended; enable filtering for true filt-skip behavior.');
end

% -------- Z-score per channel (retain sd0 for optional global epsilon) --------
X = X - mean(X,2);
sd0 = std(X,0,2); sd0(sd0==0)=1;
X = X ./ sd0;

% -------- Per-scale bands (store exact edges) --------
scales = cell(1,nScales);
for s = 1:nScales
    lo = 0; hi = nf;
    switch lower(o.FilterMode)
        case 'lowpass'
            lo = 0;          hi = nf/s;          % classic LP cascade (overlap by design)
        case 'highpass'
            lo = nf/(s+1);   hi = nf;
        case 'bandpass'
            if strcmpi(o.BandMode,'annuli')
                lo = nf/(s+1); hi = nf/s;        % adjacent non-overlapping rings
            else
                lo = max(0,o.Band(1)); hi = min(nf,o.Band(2));
            end
        case 'bandstop'
            if strcmpi(o.BandMode,'annuli')
                lo = nf/(s+1); hi = nf/s;        % stop that ring
            else
                lo = max(0,o.Band(1)); hi = min(nf,o.Band(2));
            end
        case 'none'
            lo = 0; hi = nf;
        otherwise
            error('Unknown FilterMode: %s', o.FilterMode);
    end
    lo = max(0, min(lo, nf));
    hi = max(0, min(hi, nf));
    if ~(hi > lo)
        warning('Scale %d yields invalid band [%.6g %.6g] Hz; setting NaNs.', s, lo, hi);
        scales{s} = [NaN NaN]; continue
    end
    scales{s} = [lo, hi];
    if o.Verbose
        if strcmpi(o.FilterMode,'bandstop')
            fprintf('Scale %2d: STOP [%.3f %.3f] Hz\n', s, lo, hi);
        else
            fprintf('Scale %2d: PASS [%.3f %.3f] Hz\n', s, lo, hi);
        end
    end
end

% -------- Outputs --------
mse    = nan(nChan, nScales);
cg_len = zeros(1, nScales);

% -------- Progress setup --------
useProgress = o.Progress;
useWB = useProgress && strcmpi(o.ProgressStyle,'waitbar') && usejava('desktop');
hWB = [];
if useWB
    try, hWB = waitbar(0,'Multiscale Entropy: initializing...','Name','compute_mse'); catch, hWB=[]; end
end
if useProgress && ~useWB
    fprintf('MSE: %d scales | mode=%s | coarse=%s | parallel=%s\n', nScales, o.FilterMode, ct, o.Parallel);
end
progressCount = 0; updateEvery = 1;
fmtMsg = @(c) sprintf('MSE scales: %d/%d (%.0f%%) | elapsed %s', c, nScales, 100*c/nScales, duration(0,0,toc(tStart),'format','mm:ss'));

% DataQueue for parfor progress (callback is on client only)
dq = [];
if useProgress && strcmpi(o.Parallel,'scale') && (exist('parallel.pool.DataQueue','class')==8)
    dq = parallel.pool.DataQueue;
    afterEach(dq, @(~) localUpdateProgress());
end
    function localUpdateProgress()
        progressCount = progressCount + 1;
        if mod(progressCount,updateEvery)==0
            if useWB && ~isempty(hWB) && isvalid(hWB)
                try, waitbar(progressCount/nScales, hWB, fmtMsg(progressCount)); end
            else
                fprintf('%s\n', fmtMsg(progressCount));
            end
        end
    end

% -------- Execution --------
switch lower(o.Parallel)
    case 'scale'
        MSE_s_local  = cell(1,nScales);
        cg_len_local = zeros(1,nScales);
        parfor s = 1:nScales
            lo = scales{s}(1); hi = scales{s}(2);
            Y = X;
            if filtData && ~(lo==0 && hi==nf)
                Y = mat_zpfilter(Y, lo, hi, nf, o.FiltOrder, o.FilterMode);
            end
            [CG, nBins] = coarsegrain_mat(Y, s, ct, o.CoarseMethod);
            cg_len_s = nBins;
            mse_s = nan(nChan,1);
            if nBins >= o.MinSamplesPerBin && nBins > m+1
                for ch = 1:nChan
                    cg = CG(ch,:);
                    if o.RescaleBounds
                        eps_scale = r * std(cg);
                        if eps_scale > 0, mse_s(ch) = compute_SampEn(cg, m, eps_scale, tau); else, mse_s(ch)=NaN; end
                    else
                        eps_global = r * sd0(ch);
                        mse_s(ch) = compute_SampEn(cg, m, eps_global, tau);
                    end
                end
            end
            MSE_s_local{s}  = mse_s;
            cg_len_local(s) = cg_len_s;
            if ~isempty(dq), send(dq, 1); end
        end
        for s = 1:nScales, mse(:,s) = MSE_s_local{s}; end
        cg_len = cg_len_local;

    case 'channel'
        MSE_ch_local = cell(1,nChan);
        CGL_ch_local = cell(1,nChan);
        for ch = 1:nChan
            mse_ch    = nan(1,nScales);
            cg_len_ch = zeros(1,nScales);
            for s = 1:nScales
                lo = scales{s}(1); hi = scales{s}(2);
                y = X(ch,:);
                if filtData && ~(lo==0 && hi==nf)
                    y = mat_zpfilter(y, lo, hi, nf, o.FiltOrder, o.FilterMode);
                end
                [cg_row, nBins] = coarsegrain_mat(y, s, ct, o.CoarseMethod);
                cg_len_ch(s) = nBins;
                if nBins >= o.MinSamplesPerBin && nBins > m+1
                    if o.RescaleBounds
                        eps_scale = r * std(cg_row);
                        if eps_scale > 0, mse_ch(s) = compute_SampEn(cg_row, m, eps_scale, tau); else, mse_ch(s)=NaN; end
                    else
                        eps_global = r * sd0(ch);
                        mse_ch(s) = compute_SampEn(cg_row, m, eps_global, tau);
                    end
                else
                    if o.Verbose
                        fprintf('Scale %2d skipped (ch %d): nBins=%d (< Min=%d or <= m+1=%d)\n', s, ch, nBins, o.MinSamplesPerBin, m+1);
                    end
                end
            end
            MSE_ch_local{ch} = mse_ch;
            CGL_ch_local{ch} = cg_len_ch;
            if useProgress, localUpdateProgress(); end
        end
        for ch=1:nChan, mse(ch,:) = MSE_ch_local{ch}; end
        cg_len = max(cg_len, max(cat(1,CGL_ch_local{:}),[],1)); % QC summary across chans

    otherwise % 'none'
        for s = 1:nScales
            if useProgress && ~useWB, fprintf(' - scale %d/%d\n', s, nScales); end
            lo = scales{s}(1); hi = scales{s}(2);
            Y = X;
            if filtData && ~(lo==0 && hi==nf)
                Y = mat_zpfilter(Y, lo, hi, nf, o.FiltOrder, o.FilterMode);
            end
            [CG, nBins] = coarsegrain_mat(Y, s, ct, o.CoarseMethod);
            cg_len(s) = nBins;
            if nBins >= o.MinSamplesPerBin && nBins > m+1
                for ch = 1:nChan
                    cg = CG(ch,:);
                    if o.RescaleBounds
                        eps_scale = r * std(cg);
                        if eps_scale > 0, mse(ch,s) = compute_SampEn(cg, m, eps_scale, tau); else, mse(ch,s)=NaN; end
                    else
                        eps_global = r * sd0(ch);
                        mse(ch,s) = compute_SampEn(cg, m, eps_global, tau);
                    end
                end
            else
                if o.Verbose
                    fprintf('Scale %2d skipped: nBins=%d (< Min=%d or <= m+1=%d)\n', s, nBins, o.MinSamplesPerBin, m+1);
                end
            end
            if useProgress, localUpdateProgress(); end
        end
end

% -------- Remove NaN-only scales consistently --------
nanScale = all(isnan(mse),1);
if any(nanScale)
    mse(:,nanScale) = [];
    scales(nanScale) = [];
    cg_len(nanScale) = [];
end

% -------- Close waitbar --------
if ~isempty(hWB) && isvalid(hWB), try, close(hWB); end, end

% -------- Info out --------
info = struct('FilterMode',o.FilterMode,'Band',o.Band,'RescaleBounds',o.RescaleBounds, ...
              'MinSamplesPerBin',o.MinSamplesPerBin,'FiltOrder',o.FiltOrder, ...
              'Parallel',o.Parallel,'Verbose',o.Verbose,'fs',fs, ...
              'cg_len', cg_len, 'elapsed', toc(tStart));

end

% ======================= LOCAL SUBFUNCTIONS =======================

function Y = mat_zpfilter(Y, lo, hi, nf, filtOrder, filterMode)
% Zero-phase IIR Butterworth (SOS if available) + reflection padding.
if strcmpi(filterMode,'none'), return; end
switch lower(filterMode)
    case 'lowpass',   Wn = hi/nf;       ftype='low';
    case 'highpass',  Wn = lo/nf;       ftype='high';
    case 'bandpass',  Wn = [lo hi]/nf;  ftype='bandpass';
    case 'bandstop',  Wn = [lo hi]/nf;  ftype='stop';
    otherwise, return;
end
if any(Wn<=0) || any(Wn>=1) || (numel(Wn)==2 && Wn(1)>=Wn(2)), return; end
useSOS = exist('sosfiltfilt','file')==2;
if useSOS
    [z,p,k] = butter(filtOrder, Wn, ftype); sos = zp2sos(z,p,k);
else
    [b,a] = butter(filtOrder, Wn, ftype);
end
T = size(Y,2);
pad = min(max(60, 3*filtOrder*3), floor((T-1)/2));
if pad>0
    Ypad = [fliplr(Y(:,1:pad)) , Y , fliplr(Y(:,end-pad+1:end))];
    if useSOS, Ypad = sosfiltfilt(sos, Ypad.').'; else, Ypad = filtfilt(b,a, Ypad.').'; end
    Y = Ypad(:, pad+1:end-pad);
else
    if useSOS, Y = sosfiltfilt(sos, Y.').'; else, Y = filtfilt(b,a, Y.').'; end
end
end

function [CG, nBins] = coarsegrain_mat(Y, s, ct, coarseMethod)
% Coarse-grain along time by factor s. Y: [nChan x nSamples] or [1 x nSamples].
if nargin < 4 || isempty(coarseMethod), coarseMethod = 'stat'; end
if isrow(Y), Y = Y(:)'; end
nChan = size(Y,1);
nFull = floor(size(Y,2)/s)*s;
if nFull == 0, CG = nan(nChan,0); nBins = 0; return; end

if strcmpi(coarseMethod,'skip')
    % Kosciessa filt-skip: take every s-th sample after filtering
    CG = Y(:, 1:s:nFull);
    nBins = size(CG,2);
    return
end

% Azami-style statistical coarse-grain
nBins = nFull / s;
Z = reshape(Y(:,1:nFull), nChan, s, nBins);
switch ct
    case 'mean'
        CG = squeeze(mean(Z, 2, 'omitnan'));
    case 'sd'
        CG = squeeze(std(Z, 0, 2, 'omitnan'));
    case 'variance'
        CG = squeeze(var(Z, 0, 2, 'omitnan'));
    otherwise
        error('Unknown coarseType token after parsing: %s', ct);
end
if isvector(CG), CG = CG(:)'; end
end
