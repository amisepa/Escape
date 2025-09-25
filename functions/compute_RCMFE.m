function [RCMFE, scales] = compute_RCMFE(data, varargin)
% compute_RCMFE  Composite Multiscale Fuzzy Entropy (Azami & Escudero style) via compute_FuzzEn only.
%
%   [RCMFE, scales] = compute_RCMFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                                   'n', 2, 'coarsing','std', 'num_scales', 15, ...
%                                   'MinSamplesPerBin', 4, ...
%                                   'Parallel', true, 'Progress', true)
%
% Inputs (name-value):
%   data              : EEGLAB EEG struct (EEG.data [n_ch x n_samp]) OR numeric [n_ch x n_samp]
%   'm'               : embedding dimension (default = 2)
%   'r'               : similarity bound (default = 0.15; data are z-scored per channel)
%   'tau'             : time lag (default = 1)
%   'n'               : fuzzy exponent for exponential kernel (default = 2)
%   'coarsing'        : 'median' | 'mean' | 'trimmed mean' (20%) | 'std' [default] | 'var'
%   'num_scales'      : requested number of scales (default = 15)
%   'MinSamplesPerBin': HARD guard: >= this many coarse bins to compute a scale (default = 4)
%   'Parallel'        : parfor over channels (default = true)
%   'Progress'        : header / progress reporting (default = true)
%
% Output:
%   RCMFE  : [n_channels x S] composite multiscale fuzzy entropy per channel & scale
%   scales : 1:S (simple numeric scale index; no filtering)
%
% Notes / Integrity:
%   • Signals are z-scored per channel (mean=0, std=1), so r≈0.15, n=2 are standard.
%   • Scale 1 = FuzzyEn of the original (z-scored) signal (no coarse-graining).
%   • For s>1 we coarse-grain at each offset (ii=1..s) and compute FuzzEn per coarse series,
%     then **average across offsets** (Composite MFE). This uses ONLY compute_FuzzEn().
%   • The true Refined Composite (RCMFE) requires access to pm/pm1; since compute_FuzzEn does
%     not expose those, this implementation computes the **composite** variant (CMFE).
%
% References:
%   Azami, H., Fernández, A., & Escudero, J. (2017). Med Biol Eng Comput, 55(11), 2037–2052.
%   Azami, H., & Escudero, J. (2018). Entropy, 20(2), 138.
%   Chen, W., Wang, Z., Xie, H., & Yu, W. (2007). IEEE TNSRE, 15(2), 266–272.
%   Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Phys Rev Lett, 89(6), 068102.
%   Costa, M. D., & Goldberger, A. L. (2015). Entropy, 17(3), 1197–1203.
%   Grandy, T. H., et al. (2016). Sci Rep, 6:23073.
%   Kosciessa, J. Q., et al. (2020). NeuroImage, 206:116316.  % scale treatment context; no filtering here
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','std',       @(s) any(strcmpi(s,{'median','mean','trimmed mean','trimmed','tmean','trim20','std','sd','standard deviation','var','variance'})));
p.addParameter('num_scales', 15,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
nScales_req  = p.Results.num_scales;
minBinsHard  = p.Results.MinSamplesPerBin;
parallelMode = p.Results.Parallel;
showProg     = p.Results.Progress;

% ---------- Coerce data ----------
if isstruct(data), X = double(data.data); else, X = double(data); end
if size(X,1) > size(X,2), X = X.'; end
[nch, nSamp] = size(X);

% ---------- Cap S by minimum coarse bins ----------
S    = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1, minBinsHard));
if S > maxS
    warning('compute_RCMFE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', S, maxS, minBinsHard);
    S = maxS;
end
if S < 1, S = 1; end
scales = 1:S;

% ---------- Z-score per channel (no NaN filling) ----------
Xz = X;
for c = 1:nch
    x  = X(c,:);
    mu = mean(x,'omitnan');
    sd = std(x,0,'omitnan');
    if ~isfinite(sd) || sd==0
        Xz(c,:) = 0;
    else
        Xz(c,:) = (x - mu) ./ sd;
    end
end

% ---------- Output + header ----------
RCMFE = nan(nch, S);
if showProg
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('RCMFE (composite, via compute\\_FuzzEn): %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=on\n', ...
            nch, m, tau, r, n_exp, upperLabel(coarseType), S);
    else
        fprintf('RCMFE (composite, via compute\\_FuzzEn): %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=off\n', ...
            nch, m, tau, r, n_exp, upperLabel(coarseType), S);
    end
end

% ---------- Progress UI ----------
useWB = ~parallelMode && usejava('desktop') && showProg;
hWB = [];
if useWB
    try hWB = waitbar(0,'Computing Composite MFE...','Name','compute_RCMFE'); catch, hWB = []; end
end

% For parfor: aggregate progress cleanly (no worker spam)
useDQ = parallelMode && ~isempty(ver('parallel')) && showProg;
if useDQ
    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @notifyProgress);
end

% ---------- Compute ----------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        RCMFE(ch,:) = cmfe_one_channel(Xz(ch,:), m, r, n_exp, tau, coarseType, S, minBinsHard);
        if useDQ, send(dq, 1); end
    end
else
    for ch = 1:nch
        RCMFE(ch,:) = cmfe_one_channel(Xz(ch,:), m, r, n_exp, tau, coarseType, S, minBinsHard);
        if ~isempty(hWB) && isvalid(hWB)
            try waitbar(ch/nch, hWB, sprintf('Computing Composite MFE... (%d/%d)', ch, nch)); catch, end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch, end, end
end

% ---------- nested: progress printer ----------
    function notifyProgress(~)
        nDone = nDone + 1;
        % print at start, every ~5%, and at end
        step = max(1, round(0.05*nch));
        if nDone==1 || nDone==nch || mod(nDone, step)==0
            fprintf('  progress: ch %d/%d\n', nDone, nch);
        end
    end
end % main

% ========================================================================
function v = cmfe_one_channel(sig, m, r, n_exp, tau, coarseType, S, minBinsHard)
% Composite MFE using ONLY compute_FuzzEn — mean across offsets at each scale.
v = nan(1, S);

for s = 1:S
    if s == 1
        v(1) = fuzz_one(sig, m, r, n_exp, tau);   % original series
        continue
    end

    L = floor(numel(sig)/s)*s;
    nBins = L / s;
    if nBins < max(minBinsHard, m+1), continue; end

    fe_vals = nan(1, s);
    for off = 1:s
        xoff = sig(off:end);
        Loff = floor(numel(xoff)/s)*s;
        nBins_off = Loff / s;
        if nBins_off < max(minBinsHard, m+1), continue; end

        Y = reshape(xoff(1:Loff), s, []);       % [s x nBins_off]
        cg = coarsegrain(Y, coarseType);        % 1 x nBins_off
        fe_vals(off) = fuzz_one(cg, m, r, n_exp, tau);
    end

    if any(isfinite(fe_vals))
        v(s) = mean(fe_vals(isfinite(fe_vals)));
    else
        v(s) = NaN;
    end
end
end

% ========================================================================
function cg = coarsegrain(Y, ct)
% Column-wise coarse-graining: Y is [s x nBins], output 1 x nBins
switch lower(strtrim(ct))
    case 'mean'
        cg = mean(Y, 1, 'omitnan');
    case 'median'
        cg = median(Y, 1, 'omitnan');
    case {'trimmed mean','trimmed','tmean','trim20'}
        pct = 20;   % percent total (10% per tail)
        if exist('trimmean','file')==2 && all(isfinite(Y(:)))
            cg = trimmean(Y, pct, 1);  % across rows (dim=1)
        else
            % manual trim if trimmean missing or NaNs present
            nSeg = size(Y,2); cg = nan(1,nSeg);
            kfrac = pct/200;  % 10% per tail
            for j = 1:nSeg
                col = Y(:,j); col = col(isfinite(col));
                if isempty(col), cg(j)=NaN; continue; end
                k = floor(kfrac*numel(col));
                if 2*k >= numel(col), cg(j)=NaN; else, col = sort(col); cg(j)=mean(col(k+1:end-k)); end
            end
        end
    case {'sd','std','standard deviation'}
        cg = std(Y, 0, 1, 'omitnan');
    case {'variance','var'}
        cg = var(Y, 0, 1, 'omitnan');
    otherwise
        cg = std(Y, 0, 1, 'omitnan');  % default
end
cg = cg(:).';
end

% ========================================================================
function val = fuzz_one(xrow, m, r, n_exp, tau)
% Call your compute_FuzzEn on a single row; Parallel/Progress OFF.
if ~isrow(xrow), xrow = xrow(:).'; end
val = compute_FuzzEn(double(xrow), ...
    'm', m, 'n', n_exp, 'tau', tau, 'r', r, ...
    'Kernel','exponential', 'BlockSize', 256, ...
    'Parallel', false, 'Progress', false);
if ~isscalar(val), val = val(1); end
val = double(val);
end

% ========================================================================
function s = upperLabel(coarseType)
cl = lower(strtrim(coarseType));
if any(strcmp(cl,{'sd','std','standard deviation'}))
    s = 'STD';
elseif any(strcmp(cl,{'var','variance'}))
    s = 'VAR';
elseif strcmp(cl,'mean')
    s = 'MEAN';
elseif strcmp(cl,'median')
    s = 'MEDIAN';
elseif any(strcmp(cl,{'trimmed mean','trimmed','tmean','trim20'}))
    s = 'TRIM20';
else
    s = upper(coarseType);
end
end
