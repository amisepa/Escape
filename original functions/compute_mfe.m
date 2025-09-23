function [MFE, scales] = compute_MFE(data, varargin)
% compute_MFE  Multiscale Fuzzy Entropy (Costa-style coarse-graining) across channels.
%
%   [MFE, scales] = compute_MFE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                               'n', 2, 'coarsing','std', 'num_scales', 15, ...
%                               'MinSamplesPerBin', 4, 'StableMinBins', 100, ...
%                               'Parallel', true, 'Progress', true)
%
% Inputs:
%   data              : EEGLAB EEG struct (EEG.data [n_ch x n_samp]) OR numeric [n_ch x n_samp]
%   'm'               : embedding dimension (default = 2)
%   'r'               : similarity bound (default = 0.15; data z-scored per channel)
%   'tau'             : time lag (default = 1)
%   'n'               : fuzzy exponent (default = 2)
%   'coarsing'        : 'median' | 'mean' | 'trimmed mean' (20%) | 'std' [default] | 'var'
%   'num_scales'      : requested number of scales (default = 15)
%   'MinSamplesPerBin': global guard; minimum coarse bins at every scale to keep S (default = 4)
%   'StableMinBins'   : per-scale stability floor for FE on coarse series (default = 100)
%   'Parallel'        : parfor over channels (default = true)
%   'Progress'        : print progress / waitbar (default = true)
%
% Output:
%   MFE    : [n_channels x S] multiscale fuzzy entropy per channel & scale
%   scales : 1:S (simple numeric scale index; no filtering applied)
%
% Notes:
%   • Signals are z-scored per channel (mean=0, std=1), so r≈0.15 and n=2 are standard.
%   • No bandpass/spectral filtering here.
%   • Scale 1 = FuzzyEn of the original (z-scored) signal; coarse-graining for s≥2 only.
%   • compute_fe is called with internal parallel/progress OFF (fallbacks included).
%
% References:
%   Chen, W., Wang, Z., Xie, H., & Yu, W. (2007). IEEE TNSRE, 15(2), 266–272. (FuzzyEn)
%   Azami, H., Fernández, A., & Escudero, J. (2017). Med Biol Eng Comput, 55(11), 2037–2052. (RCMFEσ/μ)
%   Azami, H., & Escudero, J. (2018). Entropy, 20(2), 138. (coarse-graining & downsampling)
%   Costa, M., Goldberger, A. L., & Peng, C.-K. (2002). Phys Rev Lett, 89(6), 068102.
%   Costa, M. D., & Goldberger, A. L. (2015). Entropy, 17(3), 1197–1203. (GMSE/volatility)
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,              @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,               @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('n', 2,                 @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('coarsing','std',       @(s) any(strcmpi(s,{'median','mean','trimmed mean','trimmed','tmean','trim20','std','sd','standard deviation','var','variance'})));
p.addParameter('num_scales', 15,       @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('MinSamplesPerBin', 4,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('StableMinBins', 100,   @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,       @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,       @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
tau          = p.Results.tau;
n_exp        = p.Results.n;
coarseType   = p.Results.coarsing;
nScales_req  = p.Results.num_scales;
minBinsAll   = p.Results.MinSamplesPerBin;
minBinsStable= p.Results.StableMinBins;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;

% ---------------- Get numeric data ----------------
if isstruct(data)
    X = double(data.data);
else
    X = double(data);
end
if size(X,1) > size(X,2), X = X.'; end
[nch, nSamp] = size(X);

% ---------------- Cap S so every scale has >= minBinsAll coarse bins -----
S    = max(1, floor(nScales_req));
maxS = floor(nSamp / max(1,minBinsAll));
if S > maxS
    warning('compute_MFE:ReducingScales', ...
        'Reducing num_scales from %d to %d to keep >=%d coarse bins at every scale.', S, maxS, minBinsAll);
    S = maxS;
end
if S < 1, S = 1; end
scales = 1:S;

% ---------------- Fill missing & z-score per channel ---------------------
for ch = 1:nch
    xi = X(ch,:);
    if any(~isfinite(xi))
        try
            xi = fillmissing(xi,'linear','EndValues','nearest');
        catch
            % fallback: forward/back fill, then average neighbors
            isn = ~isfinite(xi);
            if any(isn)
                idx = find(~isn,1,'first'); if ~isempty(idx), xi(1:idx-1) = xi(idx); end
                idx = find(~isn,1,'last');  if ~isempty(idx), xi(idx+1:end) = xi(idx); end
                prev = [xi(1),     xi(1:end-1)];
                next = [xi(2:end), xi(end)];
                bad  = isn & isfinite(prev) & isfinite(next);
                xi(bad) = 0.5*(prev(bad)+next(bad));
            end
        end
        X(ch,:) = xi;
    end
end
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

% ---------------- Outputs & progress header ------------------------------
MFE = nan(nch, S);
if showProgress
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('MFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=on\n', ...
            nch, m, tau, r, n_exp, upperLabel(coarseType), S);
    else
        fprintf('MFE: %d ch | m=%g, tau=%g, r=%g, n=%g | coarse=%s | S=%d | parallel=off\n', ...
            nch, m, tau, r, n_exp, upperLabel(coarseType), S);
    end
end

% Optional waitbar (serial mode only)
useWB = showProgress && ~parallelMode && usejava('desktop');
hWB = [];
if useWB
    try hWB = waitbar(0,'Computing Multiscale Fuzzy Entropy...','Name','compute_MFE'); catch, hWB = []; end
end

% ---------------- Compute per channel ------------------------------------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        MFE(ch,:) = mfe_one_channel(Xz(ch,:), m, r, n_exp, tau, coarseType, S, minBinsAll, minBinsStable, showProgress, ch, nch);
    end
else
    for ch = 1:nch
        MFE(ch,:) = mfe_one_channel(Xz(ch,:), m, r, n_exp, tau, coarseType, S, minBinsAll, minBinsStable, showProgress, ch, nch);
        if ~isempty(hWB) && isvalid(hWB)
            try waitbar(ch/nch, hWB, sprintf('Computing MFE... (%d/%d)', ch, nch)); catch, end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch, end, end
end
end % ---- main function end ----

% ========================================================================
function v = mfe_one_channel(sig, m, r, n_exp, tau, coarseType, S, minBinsAll, minBinsStable, showProgress, ch, nch)
% Compute MFE across scales (no filtering)

v = nan(1, S);
for s = 1:S
    if s == 1
        % Scale 1 = FuzzyEn of original (z-scored) signal
        v(1) = fe_safe(sig, m, r, n_exp, tau);
        continue
    end

    L = floor(numel(sig)/s)*s;
    nBins = L / s;
    minNeeded = max([minBinsAll, m+1, minBinsStable]);

    if nBins < minNeeded
        if showProgress
            fprintf('  [drop] ch %d: scale %d nBins=%d (<%d)\n', ch, s, nBins, minNeeded);
        end
        continue
    end

    y = reshape(sig(1:L), s, []);   % [s x nBins]
    switch lower(coarseType)
        case 'mean'
            cg = mean(y, 1, 'omitnan');
        case 'median'
            cg = median(y, 1, 'omitnan');
        case {'trimmed mean','trimmed','tmean','trim20'}
            pct = 0.20;
            if exist('trimmean','file')==2 && all(isfinite(y(:)))
                cg = trimmean(y, pct*100);
            else
                nSeg = size(y,2); cg = nan(1,nSeg);
                for j = 1:nSeg
                    col = y(:,j); col = col(isfinite(col));
                    if isempty(col), cg(j)=NaN; continue; end
                    k = floor((pct/2)*numel(col)); % 10% per tail
                    if 2*k >= numel(col), cg(j)=NaN; else, col = sort(col); cg(j)=mean(col(k+1:end-k)); end
                end
            end
        case {'sd','std','standard deviation'}
            cg = std(y, 0, 1, 'omitnan');
        case {'variance','var'}
            cg = var(y, 0, 1, 'omitnan');
        otherwise
            cg = std(y, 0, 1, 'omitnan');
    end
    cg = cg(:).';
    v(s) = fe_safe(cg, m, r, n_exp, tau);
end

if showProgress
    fprintf('  ch %3d/%3d: done\n', ch, nch);
end
end

% ========================================================================
function val = fe_safe(x, m, r, n_exp, tau)
% Try compute_fe / compute_FE / FuzzEn-style functions; no internal parallel/progress flags
try
    val = compute_fe(x, m, r, n_exp, tau);
catch
    try
        val = compute_FE(x, m, r, n_exp, tau);
    catch
        try
            % if a single-channel FuzzyEn function exists (e.g., fuzz), try it
            val = fuzz(x, m, r, n_exp, tau);
        catch
            try
                % fall back to a blockwise fuzzy entropy if available
                val = fuzz_blockwise(x, m, r, n_exp, tau, 'exponential', 256);
            catch
                error('compute_MFE: Fuzzy entropy function not found (compute_fe / compute_FE / fuzz[_blockwise]).');
            end
        end
    end
end
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
