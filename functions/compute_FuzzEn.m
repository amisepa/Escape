function FuzzEn = compute_FuzzEn(data, varargin)
% compute_FuzzEn  Computes Fuzzy Entropy (FuzzEn) across multichannel data.
%
%   FuzzEn = compute_FuzzEn(data, 'm', 2, 'n', 2, 'tau', 1, 'r', .15, ...
%                           'Kernel','exponential', 'BlockSize', 256, 
%                           'Parallel', false, 'Progress', false)
%
% Inputs:
%   data        : EEG data matrix [n_channels x n_samples] (numeric)
%   'm'         : embedding dimension (default = 2)
%   'n'         : fuzziness exponent for exponential kernel (default = 2)
%   'tau'       : time lag (default = 1)
%   'r'         : similarity bound (default = .15)
%   'Kernel'    : 'exponential' [default] | 'gaussian'
%   'BlockSize' : pair-block size to bound memory (default = 256; >=128)
%   'Parallel'  : logical true/false to enable parfor over channels (default = false)
%   'Progress'  : logical true/false to show progress (default = false)
%
% Output:
%   FuzzEn      : [n_channels x 1] Fuzzy Entropy per channel
%
% Notes:
%   • Data is z-score normalized per channel so std = 1
%   • Uses bounded-memory, blockwise pair counting (no full M×M broadcasts).
%
% References:
%   Chen, W., Wang, Z., Xie, H., & Yu, W. (2007). Characterization of 
%       surface EMG signal based on fuzzy entropy. IEEE Transactions on 
%       neural systems and rehabilitation engineering, 15(2), 266-272.
%   Azami, H., Fernández, A., & Escudero, J. (2017). Refined multiscale 
%       fuzzy entropy based on standard deviation for biomedical signal 
%       analysis. Medical & biological engineering & computing, 55(11), 2037-2052.
% 
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,               @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('n', 2,               @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1,             @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', .15,             @(x) isnumeric(x) && x > 0 && x < 1);
p.addParameter('Kernel','exponential',@(s) any(strcmpi(s,{'exponential','gaussian'})));
p.addParameter('BlockSize', 256,    @(x) isnumeric(x) && isscalar(x) && x >= 128 && x <= 2048);
p.addParameter('Parallel', true,    @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true,    @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
n_exp        = p.Results.n;
tau          = p.Results.tau;
r            = p.Results.r;
kernelType   = lower(p.Results.Kernel);
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;
blk          = p.Results.BlockSize;

if size(data,1) > size(data,2)
    data = data.';  % [n_channels x n_samples]
end
[nchan, ~] = size(data);

% Z-score normalization across time per channel (Azami et al., 2020)
% data_z = normalize(data, 2);    % one line but introduced in 2018a
data_z = data;
for c = 1:nchan
    x  = data(c,:);
    sd = std(x,0,'omitnan');
    if ~isfinite(sd) || sd == 0
        data_z(c,:) = 0;
    else
        mu = mean(x,'omitnan');
        data_z(c,:) = (x - mu)./sd;
    end
end



FuzzEn = nan(nchan, 1);

%  Progress headers 
if showProgress
    if parallelMode
        fprintf('FuzzEn: %d channel(s) | m=%g, tau=%g, r=%g, n=%g, | kernel=%s | parallel=on (text only)\n', ...
            nchan, m, tau, r, n_exp, kernelType);
        fprintf('Progress:\n');
    else
        fprintf('FuzzEn: %d channel(s) | m=%g, tau=%g, r=%g, n=%g,  | kernel=%s | parallel=off (text + waitbar)\n', ...
            nchan, m, tau, r, n_exp, kernelType);
    end
end


%  Compute per channel 
if parallelMode && ~isempty(ver('parallel'))
    parfor iChan = 1:nchan
        FuzzEn(iChan) = fuzz_blockwise(data_z(iChan,:), m, r, n_exp, tau, kernelType, blk);
        % FuzzEn(iChan) = fuzz(data_z(iChan,:), m, r, n_exp, tau, kernelType);
        % if showProgress && (mod(iChan, max(1, floor(nchan/20)))==0 || iChan==nchan)
        if showProgress
            fprintf('  ch %3d/%3d: %.6f\n', iChan, nchan, FuzzEn(iChan));
        end
    end
else
    useWB = showProgress && usejava('desktop');
    hWB = [];
    if useWB
        try hWB = waitbar(0,'Computing Fuzzy Entropy...','Name','compute_FuzzEn'); catch, hWB = []; end
    end
    for iChan = 1:nchan
        FuzzEn(iChan) = fuzz_blockwise(data_z(iChan,:), m, r, n_exp, tau, kernelType, blk);
        % FuzzEn(iChan) = fuzz(data_z(iChan,:), m, r, n_exp, tau, kernelType);
        if showProgress
            fprintf('  ch %3d/%3d: %.6f\n', iChan, nchan, FuzzEn(iChan));
            if ~isempty(hWB) && isvalid(hWB)
                try waitbar(iChan/nchan, hWB, sprintf('Computing Fuzzy Entropy... (%d/%d)', iChan, nchan)); catch; end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch; end; end
end
end

% =========================================================================
function fe = fuzz_blockwise(signal, m, r, n, tau, kernelType, blk)
% function fe = fuzz(signal, m, r, n, tau, kernelType)
% Bounded-memory FuzzyEn with Chebyshev metric.

x = signal(isfinite(signal));
if tau > 1, x = x(1:tau:end); end
N = numel(x);
if N <= m + 1, fe = NaN; return; end

% Embeddings for k=m and k=m+1
[Xm,  Lm]  = embed_tau(x, m,   tau);
[Xm1, Lm1] = embed_tau(x, m+1, tau);
if Lm < 2 || Lm1 < 2, fe = NaN; return; end

% Mean fuzzy similarity over UNIQUE pairs
pm  = fuzzy_mean_similarity_block(Xm,  r, n, kernelType, blk);
pm1 = fuzzy_mean_similarity_block(Xm1, r, n, kernelType, blk);
% pm  = fuzzy_mean_similarity_fast(Xm,  r, n, kernelType);
% pm1 = fuzzy_mean_similarity_fast(Xm1, r, n, kernelType);

if pm>0 && pm1>0 && isfinite(pm) && isfinite(pm1)
    fe = log(pm / pm1);
else
    fe = NaN;
end
end

function mu_mean = fuzzy_mean_similarity_block(X, r, n, kernelType, blk)
% Mean fuzzy similarity over unique pairs; diagonal uses upper triangle only.

nVec = size(X,1);
if nVec < 2, mu_mean = NaN; return; end

sum_mu = 0; num_p = 0;

for i1 = 1:blk:nVec
    i2 = min(i1+blk-1, nVec);
    Xi = X(i1:i2,:);  bi = size(Xi,1);

    % (A) Diagonal part (upper triangle)
    if bi > 1
        for d1 = 1:blk:bi
            d2 = min(d1+blk-1, bi);
            Xd = Xi(d1:d2,:); bd = size(Xd,1);
            if bd > 1
                Dmax = zeros(bd,bd);
                for dim = 1:size(X,2)
                    Dij = abs(Xd(:,dim) - Xd(:,dim).');   % bi-di distances in this sub-block
                    if dim==1, Dmax = Dij; else, Dmax = max(Dmax, Dij); end
                end
                ut = triu(true(bd,bd),1);
                dvec = Dmax(ut);
                sum_mu = sum_mu + sum( fuzzy_kernel(dvec, r, n, kernelType) );
                num_p  = num_p  + nnz(ut);
            end
        end
    end

    % (B) Off-diagonals (full bi×bj)
    for j1 = (i2+1):blk:nVec
        j2 = min(j1+blk-1, nVec);
        Xj = X(j1:j2,:);  bj = size(Xj,1);
        Dmax = zeros(bi,bj);
        for dim = 1:size(X,2)
            Dij = abs(Xi(:,dim) - Xj(:,dim).');
            if dim==1, Dmax = Dij; else, Dmax = max(Dmax, Dij); end
        end
        sum_mu = sum_mu + sum( fuzzy_kernel(Dmax(:), r, n, kernelType) );
        num_p  = num_p  + bi*bj;
    end
end

mu_mean = tern(num_p==0, NaN, sum_mu/num_p);
end

function y = fuzzy_kernel(d, r, n, kernelType)
switch kernelType
    case 'exponential', y = exp(-(d.^n) / r);
    case 'gaussian',    y = exp(-(d.^2) / (2*r^2));
    otherwise,          y = exp(-(d.^n) / r);
end
end

function [X, L] = embed_tau(signal, m, tau)
signal = signal(:).';
N = numel(signal);
L = N - (m-1)*tau;
if L <= 0, X = zeros(0,m); return; end
idx  = (0:(m-1))*tau;
rows = (1:L).';
X = signal(rows + idx);
end

function y = tern(cond, a, b), if cond, y=a; else, y=b; end, end
