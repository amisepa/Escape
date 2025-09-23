function [MSE, scales] = compute_MSE(data, varargin)
% compute_MSE  Multiscale Entropy (Costa-style) across multichannel data.
%
%   [MSE, scales] = compute_MSE(data, 'm', 2, 'r', 0.15, 'tau', 1, ...
%                               'CoarseType','SD', 'nScales', 15, ...
%                               'Parallel', true, 'Progress', true)
%
% Inputs:
%   data        : EEGLAB EEG struct (uses EEG.data [n_channels x n_samples])
%                 OR numeric matrix [n_channels x n_samples] (preferred orientation)
%   'm'         : embedding dimension (default = 2)
%   'r'         : similarity bound (default = 0.15; data are z-scored per channel)
%   'tau'       : time lag (default = 1)
%   'coarsing'  : 'median' (default), 'mean', 'trimmed mean', 'std' [default], 'var'
%   'nScales'   : number of scales (default = 15)
%   'Parallel'  : logical true/false, parfor over channels (default = true)
%   'Progress'  : logical true/false, print progress (default = true)
%
% Output:
%   MSE         : [n_channels x nScales] multiscale entropy per channel & scale
%   scales      : 1:nScales (simple numeric scale index; no filtering applied)
%
% Notes:
%   • Signals are z-scored per channel (mean=0, std=1), so r≈0.15 is standard.
%   • No bandpass or spectral filtering is performed here.
%   • Internally calls compute_sampEn with Parallel/Progress disabled.
%
% References (kept from your doc):
%   Azami, H., & Escudero, J. (2018). Coarse-graining approaches in 
%       univariate multiscale sample and dispersion entropy. Entropy, 20(2), 138.
%   Azami, H., Fernández, A., & Escudero, J. (2017). Refined multiscale 
%       fuzzy entropy based on standard deviation for biomedical signal 
%       analysis. Medical & biological engineering & computing, 55(11), 2037-2052.
%   Costa, M., Goldberger, A. L., & Peng, C.-K. (2002).
%       Multiscale Entropy Analysis of Complex Physiologic Time Series.
%       Phys. Rev. Lett., 89(6), 068102.
%   Costa, M. D., & Goldberger, A. L. (2015). Generalized multiscale
%       entropy analysis: Application to quantifying the complex volatility
%       of human heartbeat time series. Entropy, 17(3), 1197-1203.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
% data can be EEG struct or numeric [nCh x nSamp]
p.addRequired('data', @(x) (isstruct(x) && isfield(x,'data')) || (isnumeric(x) && ndims(x)==2));
p.addParameter('m', 2,                @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('r', 0.15,             @(x) isnumeric(x) && isscalar(x) && x>0 && x<2);
p.addParameter('tau', 1,              @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('coarsing','std',      @(s) any(strcmpi(s,{'median', 'mean', 'trimmed mean', 'std','var'})));
p.addParameter('num_scales', 15,      @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Parallel', true,      @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,      @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r            = p.Results.r;
tau          = p.Results.tau;
coarseType   = p.Results.coarsing;
nScales      = p.Results.num_scales;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;

% ---------------- Get numeric data ----------------
if isstruct(data)
    X = data.data;   % [nCh x nSamples]
else
    X = data;
end
if size(X,1) > size(X,2)
    X = X.'; % enforce [nCh x nSamp]
end
[nch, ~] = size(X);

% ---------------- Z-score per channel ----------------
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

% ---------------- Outputs ----------------
MSE    = nan(nch, nScales);
scales = 1:nScales;

% ---------------- Progress header ----------------
if showProgress
    if parallelMode && ~isempty(ver('parallel'))
        fprintf('MSE: %d channel(s) | m=%g, tau=%g, r=%g | coarse=%s | nScales=%d | parallel=on\n', ...
            nch, m, tau, r, upperLabel(coarseType), nScales);
    else
        fprintf('MSE: %d channel(s) | m=%g, tau=%g, r=%g | coarse=%s | nScales=%d | parallel=off\n', ...
            nch, m, tau, r, upperLabel(coarseType), nScales);
    end
end

% ---------------- Compute per channel ----------------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nch
        MSE(ch,:) = mse_one_channel(Xz(ch,:), m, r, tau, coarseType, nScales, showProgress, ch, nch);
    end
else
    useWB = showProgress && usejava('desktop');
    hWB = [];
    if useWB
        try hWB = waitbar(0,'Computing Multiscale Entropy...','Name','compute_MSE'); catch, hWB = []; end
    end
    for ch = 1:nch
        MSE(ch,:) = mse_one_channel(Xz(ch,:), m, r, tau, coarseType, nScales, showProgress, ch, nch);
        if ~isempty(hWB) && isvalid(hWB)
            try waitbar(ch/nch, hWB, sprintf('Computing MSE... (%d/%d)', ch, nch)); catch, end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch, end, end
end

end % ---- main function end ----

% ========================================================================
function v = mse_one_channel(sig, m, r, tau, coarseType, nScales, showProgress, ch, nch)
v = nan(1, nScales);
for s = 1:nScales

    % ----- NEW: make scale 1 = plain SampEn of the (z-scored) signal -----
    if s == 1
        try
            v(1) = compute_sampEn(sig, 'm', m, 'r', r, 'tau', tau, 'Parallel', false, 'Progress', false);
        catch
            try
                v(1) = compute_sampEn(sig, m, r, tau);
            catch
                try
                    v(1) = compute_SampEn(sig, 'm', m, 'r', r, 'tau', tau, 'Parallel', false, 'Progress', false);
                catch
                    try
                        v(1) = compute_SampEn(sig, m, r, tau);
                    catch
                        v(1) = compute_se(sig, m, r, tau);
                    end
                end
            end
        end
        continue
    end
    % ---------------------------------------------------------------------

    L = floor(numel(sig)/s)*s;
    if L < max(4, (m+1)*tau)
        v(s) = NaN; continue;
    end
    y = reshape(sig(1:L), s, []);

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
                    k = floor((pct/2)*numel(col));
                    if 2*k >= numel(col), cg(j)=NaN; else, col = sort(col); cg(j)=mean(col(k+1:end-k)); end
                end
            end
        case {'sd','std','standard deviation'}
            cg = std(y, 0, 1, 'omitnan');
        case {'variance','var'}
            cg = var(y, 0, 1, 'omitnan');
        otherwise
            cg = median(y, 1, 'omitnan');  % default = median
    end
    cg = cg(:).';

    % SampEn on coarse-grained series
    try
        v(s) = compute_sampEn(cg, 'm', m, 'r', r, 'tau', tau, 'Parallel', false, 'Progress', false);
    catch
        try
            v(s) = compute_sampEn(cg, m, r, tau);
        catch
            try
                v(s) = compute_SampEn(cg, 'm', m, 'r', r, 'tau', tau, 'Parallel', false, 'Progress', false);
            catch
                try
                    v(s) = compute_SampEn(cg, m, r, tau);
                catch
                    v(s) = compute_se(cg, m, r, tau);
                end
            end
        end
    end
end



if showProgress
    fprintf('  ch %3d/%3d: done\n', ch, nch);
end
end

% ========================================================================
function s = upperLabel(coarseType)
% small label helper for header only (inline to keep file self-contained)
if strcmpi(coarseType,'SD') || strcmpi(coarseType,'Standard deviation')
    s = 'SD';
elseif strcmpi(coarseType,'Mean')
    s = 'Mean';
elseif strcmpi(coarseType,'Variance')
    s = 'Variance';
else
    s = char(coarseType);
end
end
