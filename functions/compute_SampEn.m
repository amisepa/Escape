function SampEn = compute_SampEn(data, varargin)
% compute_SampEn  Computes Sample Entropy (SampEn) across multichannel data.
%
%   SampEn = compute_SampEn(data, 'm', 2, 'tau', 1, 'r', .15, 'Parallel', false, 'Progress', false)
%
% Inputs:
%   data        : EEG data matrix [n_channels x n_samples]
%   'm'         : embedding dimension (default = 2)
%   'tau'       : time lag (default = 1)
%   'r'         : similarity bound (default = .15)
%   'Parallel'  : logical true/false to enable parfor over channels (default = false)
%   'Progress'  : logical true/false to show progress
%                 • If Parallel==true  and Progress==true → text only (parfor-safe)
%                 • If Parallel==false and Progress==true → text + waitbar (fallback to text if headless)
%
% Output:
%   SampEn      : [n_channels x 1] Sample Entropy per channel
%
% Notes:
%   • Data is z-scored per channel; r = 0.15 * std(z) per channel (Azami).
%   • Printing inside parfor may appear out of order (expected).
% 
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Ascent plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2,   @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', .15, @(x) isnumeric(x) && x > 0 && x < 1);
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
tau          = p.Results.tau;
r            = p.Results.r;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;

if size(data,1) > size(data,2)
    data = data.';  % [n_channels x n_samples]
end
[nchan, ~] = size(data);

% Z-score normalization across time per channel (Azami et al., 2020)
data_z = data;
for c = 1:nchan
    x = data(c,:);
    mu = mean(x,'omitnan');
    sd = std(x,0,'omitnan');
    if ~isfinite(sd) || sd == 0
        data_z(c,:) = 0;      % flat or all-NaN → zero; SampEn will NaN later if too short
    else
        data_z(c,:) = (x - mu)./sd;
    end
end


SampEn = nan(nchan, 1);

% Progress headers 
if showProgress
    if parallelMode
        fprintf('SampEn: %d channel(s) | m=%g, tau=%g, r =%g | parallel=on (text only)\n', nchan, m, tau, r);
    else
        fprintf('SampEn: %d channel(s) | m=%g, tau=%g, r =%g | parallel=off (text + waitbar)\n', nchan, m, tau, r);
    end
end

%  Compute per channel 
if parallelMode && ~isempty(ver('parallel'))
    % PARFOR (text-only progress if requested)
    if showProgress, fprintf('Progress:\n'); end
    parfor iChan = 1:nchan
        SampEn(iChan) = compute_SampEn_single(data_z(iChan,:), m, r, tau);
        if showProgress
            % Minimal single-line prints to reduce interleaving noise
            fprintf('  ch %3d/%3d: %.6f\n', iChan, nchan, SampEn(iChan));
        end
    end
else
    % SERIAL (text + waitbar if requested and desktop available)
    useWB = showProgress && usejava('desktop');
    hWB = [];
    if useWB
        try hWB = waitbar(0,'Computing Sample Entropy...','Name','compute_SampEn'); catch, hWB = []; end
    end
    for iChan = 1:nchan
        SampEn(iChan) = compute_SampEn_single(data_z(iChan,:), m, r, tau);
        if showProgress
            fprintf('  ch %3d/%3d: %.6f\n', iChan, nchan, SampEn(iChan));
            if ~isempty(hWB) && isvalid(hWB)
                try waitbar(iChan/nchan, hWB, sprintf('Computing Sample Entropy... (%d/%d)', iChan, nchan)); catch; end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try close(hWB); catch; end, end
end
end


% =========================================================================

function SampEn = compute_SampEn_single(signal, m, r, tau)
% Ensure finite samples only (drop NaNs/Infs), preserve order
signal = signal(isfinite(signal));
N = numel(signal);
if N <= m + 1
    SampEn = NaN; return;
end

try
    % Delay embeddings with tau
    Xm  = embed_tau(signal, m,   tau);   % [L_m  x m]
    Xm1 = embed_tau(signal, m+1, tau);   % [L_m1 x (m+1)]
    Lm  = size(Xm,  1);
    Lm1 = size(Xm1, 1);

    if Lm < 2 || Lm1 < 2
        SampEn = NaN; return;
    end

    % Chebyshev distance over unique pairs
    Dm  = pdist(Xm,  'chebychev');
    Dm1 = pdist(Xm1, 'chebychev');

    % Match probabilities (≤ r or < r; choose one consistently)
    % Use <= r by default:
    Bm  = mean(Dm  <= r);   % matches for length m
    Am1 = mean(Dm1 <= r);   % matches for length m+1

    if Bm > 0 && Am1 > 0
        SampEn = -log(Am1 / Bm);
    else
        SampEn = NaN;
    end
catch
    % Fallback: blockwise counting with same embedding & normalization
    SampEn = fallback_SampEn(signal, m, r, tau);
end
end



function SampEn = fallback_SampEn(signal, m, r, tau)
% Robust but slower; identical definition as fast path.
signal = signal(isfinite(signal));
N = numel(signal);
if N <= m + 1
    SampEn = NaN; return;
end

Xm  = embed_tau(signal, m,   tau);
Xm1 = embed_tau(signal, m+1, tau);
Lm  = size(Xm,  1);
Lm1 = size(Xm1, 1);

if Lm < 2 || Lm1 < 2
    SampEn = NaN; return;
end

% Count matches over unique pairs without pdist (blockwise to save RAM)
% Helper to count <= r Chebyshev matches
    function p = match_prob(X)
        n  = size(X,1);
        tp = n*(n-1)/2;
        if tp == 0, p = 0; return; end
        % Block size tradeoff: memory vs speed
        blk = 1024;
        hits = 0;
        for i1 = 1:blk:n-1
            j1 = min(i1+blk-1, n-1);
            Xi = X(i1:j1,:);                 % [bi x m]
            for j = (i1+1):blk:n
                j2 = min(j+blk-1, n);
                Xj = X(j:j2,:);              % [bj x m]
                % Chebyshev distance: max(abs(Xi - Xj))
                % Compute all pairwise via broadcasting:
                % diff shape -> [bi x bj x m]
                % We avoid 3D arrays: do per-dimension max incrementally
                bi = size(Xi,1); bj = size(Xj,1); mm = size(X,2);
                maxd = zeros(bi,bj);
                for d = 1:mm
                    Dij = abs(Xi(:,d) - Xj(:,d).');
                    if d == 1
                        maxd = Dij;
                    else
                        maxd = max(maxd, Dij);
                    end
                end
                hits = hits + nnz(maxd <= r);
            end
        end
        p = hits / tp;
    end

Bm  = match_prob(Xm);
Am1 = match_prob(Xm1);

if Bm > 0 && Am1 > 0
    SampEn = -log(Am1 / Bm);
else
    SampEn = NaN;
end
end


function X = embed_tau(signal, m, tau)
% Build delay-embedded vectors with step = tau (Chebyshev norm uses columns)
% Returns [n_vec x m] with rows = embedded vectors
signal = signal(:).';                       % row
N = numel(signal);
L = N - (m-1)*tau;                          % number of vectors
if L <= 0, X = zeros(0,m); return; end
idx = (0:(m-1)) * tau;
rows = (1:L).';
X = signal(rows + idx);
end