function [HD, HA, HDA, info] = compute_ExSEnt2(data, varargin)
% compute_ExSEnt  Extrema-Segmented Entropy (ExSEnt) for multichannel data.
%
%   [HD, HA, HDA, info] = compute_ExSEnt(data, 'm', 2, ...
%                                        'r', 0.20, 'lambda', 0.01, ...
%                                        'Parallel', true, 'Progress', true, ...
%                                        'Plot', false)
%
% Inputs:
%   data        : matrix [n_channels x n_samples]
%   'm'         : embedding dimension (default = 2)
%   'r'         : r scaling (r = r * std), default = 0.20
%   'lambda'    : threshold factor for extrema detection (θ = λ * IQR(diff(x))), default = 0.01
%   'Parallel'  : logical true/false to enable parfor over channels (default = true)
%   'Progress'  : logical true/false to show progress (default = true)
%                 • If Parallel==true  and Progress==true → text only (parfor-safe)
%                 • If Parallel==false and Progress==true → text + waitbar (fallback to text if headless)
%   'Plot'      : validation plot for the FIRST channel (default = false)
%
% Outputs:
%   HD, HA, HDA : [n_channels x 1] entropies (duration, amplitude, joint)
%   info        : struct with per-channel diagnostics:
%                 .M, .rD, .rA, .rDA, .rangeD, .rangeA, .segment_idx
%
% Notes:
%   • Durations (D) are in samples; amplitudes (A) are net differences.
%   • Joint SampEn uses interleaved z-scored [D A] (embedding = 2*m).
%   • Uses the SAME SampEn backend as compute_SampEn (embed_tau + Chebyshev).
%
% Example:
%   [HD,HA,HDA,info] = compute_ExSEnt(EEG, 'm', 2, 'r', 0.2, ...
%                                     'lambda', 0.01, 'Parallel', true, ...
%                                     'Progress', true, 'Plot', true);
%
% Reference:
%   Kamali, S., Baroni, F., & Varona, P. (2025). ExSEnt: Extrema-Segmented
%   Entropy Analysis of Time Series. arXiv:2509.07751
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('r', 0.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('lambda', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Plot', false, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m            = p.Results.m;
r        = p.Results.r;
lambda       = p.Results.lambda;
parallelMode = p.Results.Parallel;
showProgress = p.Results.Progress;
doPlot       = p.Results.Plot;

% ---------------- Shape ----------------
if size(data,1) > size(data,2)
    data = data.'; % [n_channels x n_samples]
end
[nchan, ~] = size(data);

HD  = nan(nchan,1);
HA  = nan(nchan,1);
HDA = nan(nchan,1);
info = struct('M',nan(nchan,1),'rD',nan(nchan,1),'rA',nan(nchan,1), ...
              'rDA',nan(nchan,1),'rangeD',nan(nchan,1),'rangeA',nan(nchan,1), ...
              'segment_idx',[]);

% ---------------- Progress header ----------------
if showProgress
    if parallelMode
        fprintf('ExSEnt: %d channel(s) | m=%g, r=%.3f, lambda=%.4f | parallel=on (text only)\n', ...
                nchan, m, r, lambda);
        fprintf('Progress:\n');
    else
        fprintf('ExSEnt: %d channel(s) | m=%g, r=%.3f, lambda=%.4f | parallel=off (text + waitbar)\n', ...
                nchan, m, r, lambda);
    end
end

% ---------------- Compute per channel ----------------
if parallelMode && ~isempty(ver('parallel'))
    parfor ch = 1:nchan
        [HD(ch),HA(ch),HDA(ch),info(ch)] = exsent_channel(data(ch,:), m, r, lambda);
        if showProgress
            fprintf('  ch %3d/%3d: M=%4d | HD=%.6f HA=%.6f HDA=%.6f\n', ...
                    ch, nchan, info(ch).M, HD(ch), HA(ch), HDA(ch));
        end
    end
else
    useWB = showProgress && usejava('desktop');
    hWB = [];
    if useWB
        try hWB = waitbar(0,'Computing ExSEnt...','Name','compute_ExSEnt'); catch, hWB = []; end
    end
    for ch = 1:nchan
        [HD(ch),HA(ch),HDA(ch),info(ch)] = exsent_channel(data(ch,:), m, r, lambda);
        if showProgress
            fprintf('  ch %3d/%3d: M=%4d | HD=%.6f HA=%.6f HDA=%.6f\n', ...
                    ch, nchan, info(ch).M, HD(ch), HA(ch), HDA(ch));
            if ~isempty(hWB) && isvalid(hWB)
                try, waitbar(ch/nchan, hWB, sprintf('Computing ExSEnt... (%d/%d)', ch, nchan)); end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); end, end
end

% ---------------- Validation plot (first channel) ----------------
if doPlot && nchan >= 1
    x = data(1,:);
    [D_vals, A_vals, segment_idx] = extract_DA(x, lambda);

    figure('Name','ExSEnt Validation','Color','w');
    subplot(2,1,1);
    plot(x,'k'); hold on;
    plot(segment_idx, x(segment_idx),'ro','MarkerFaceColor','r');
    xlabel('Samples'); ylabel('Amplitude');
    title('Signal with detected extrema');
    grid on;

    subplot(2,1,2);
    yyaxis left;  plot(D_vals,'b-o','LineWidth',1.2); ylabel('Durations (samples)');
    yyaxis right; plot(A_vals,'r-s','LineWidth',1.2); ylabel('Amplitudes (Δ)');
    xlabel('Segment index');
    title('Duration and amplitude sequences (D and A)');
    grid on;
end
end

% ========================= Channel-level helper ========================= %
function [HD,HA,HDA,info] = exsent_channel(x, m, r, lambda)
% Segment signal into durations/amplitudes
[D_vals, A_vals, segment_idx] = extract_DA(x, lambda);
M = numel(D_vals);

% Clean non-finite / zero-variance edge cases
D_vals = D_vals(isfinite(D_vals));
A_vals = A_vals(isfinite(A_vals));

if M <= m+1 || numel(D_vals) <= m+1 || numel(A_vals) <= m+1 ...
   || std(D_vals) == 0 || std(A_vals) == 0
    [HD,HA,HDA] = deal(NaN);
    info = struct('M',M,'rD',NaN,'rA',NaN,'rDA',NaN, ...
                  'rangeD',range_safe(D_vals),'rangeA',range_safe(A_vals), ...
                  'segment_idx',segment_idx);
    return;
end

% Tolerances (Azami-style on raw sequences)
rD = r * std(D_vals);
rA = r * std(A_vals);

% Use the SAME SampEn backend as compute_SampEn (embed_tau + Chebyshev)
HD = compute_SampEn_single(D_vals, m, rD, 1);
HA = compute_SampEn_single(A_vals, m, rA, 1);

% Normalize then interleave [D_norm, A_norm] → 1D series, embed with 2*m
D_norm = normalize(D_vals);
A_norm = normalize(A_vals);
joint_data = reshape([D_norm(:) A_norm(:)]', [], 1);

% Joint r after z-scoring (std≈1): use r directly
HDA = compute_SampEn_single(joint_data, 2*m, r, 1);

% Diagnostics
info = struct('M',M,'rD',rD,'rA',rA,'rDA',r, ...
              'rangeD',range_safe(D_vals),'rangeA',range_safe(A_vals), ...
              'segment_idx',segment_idx);
end

% ========================= Segmentation helper ========================== %
function [D_vals, A_vals, segment_idx] = extract_DA(signal, lambda)
% Extract segment durations and net amplitudes via robust sign-change rule.
if nargin < 2, lambda = 0.001; end
signal = signal(:)';                 % row

dx = diff(signal);
threshold = lambda * iqr(dx);

D_vals = [];
A_vals = [];
segment_idx = [];

start_idx = 1;
for i = 2:length(dx)
    if sign(dx(i)) ~= sign(dx(i-1)) && abs(dx(i) - dx(i-1)) > threshold
        duration  = i - start_idx;
        amplitude = signal(i) - signal(start_idx);

        D_vals(end+1)      = duration;   %#ok<AGROW>
        A_vals(end+1)      = amplitude;  %#ok<AGROW>
        segment_idx(end+1) = i;          %#ok<AGROW>

        start_idx = i;
    end
end
end

% ========================= SampEn backend (same as compute_SampEn) ====== %
function SampEn = compute_SampEn_single(signal, m, r, tau)
% Fast SampEn using embed_tau + Chebyshev distance; NaN-robust.
signal = signal(isfinite(signal));
N = numel(signal);
if N <= m + 1
    SampEn = NaN; return;
end

try
    Xm  = embed_tau(signal, m,   tau);   % [L_m  x m]
    Xm1 = embed_tau(signal, m+1, tau);   % [L_m1 x (m+1)]
    Lm  = size(Xm,  1);
    Lm1 = size(Xm1, 1);
    if Lm < 2 || Lm1 < 2
        SampEn = NaN; return;
    end

    Dm  = pdist(Xm,  'chebychev');
    Dm1 = pdist(Xm1, 'chebychev');

    Bm  = mean(Dm  <= r);   % matches at length m
    Am1 = mean(Dm1 <= r);   % matches at length m+1

    if Bm > 0 && Am1 > 0
        SampEn = -log(Am1 / Bm);
    else
        SampEn = NaN;
    end
catch
    SampEn = fallback_SampEn(signal, m, r, tau);
end
end

function SampEn = fallback_SampEn(signal, m, r, tau)
% Robust fallback with blockwise counting (Chebyshev).
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

    function p = match_prob(X)
        n  = size(X,1);
        tp = n*(n-1)/2;
        if tp == 0, p = 0; return; end
        blk = 1024;
        hits = 0;
        for i1 = 1:blk:n-1
            j1 = min(i1+blk-1, n-1);
            Xi = X(i1:j1,:);                 
            for j = (i1+1):blk:n
                j2 = min(j+blk-1, n);
                Xj = X(j:j2,:);              
                bi = size(Xi,1); bj = size(Xj,1); mm = size(X,2);
                maxd = zeros(bi,bj);
                for d = 1:mm
                    Dij = abs(Xi(:,d) - Xj(:,d).');
                    if d == 1, maxd = Dij; else, maxd = max(maxd, Dij); end
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
% Delay-embedded vectors with step = tau; returns [n_vec x m] (rows=templates)
signal = signal(:).';                       
N = numel(signal);
L = N - (m-1)*tau;                          
if L <= 0, X = zeros(0,m); return; end
idx  = (0:(m-1)) * tau;
rows = (1:L).';
X = signal(rows + idx);
end

% ========================= misc ========================= %
function r = range_safe(x)
x = x(isfinite(x));
if isempty(x), r = NaN; else, r = max(x) - min(x); end
end
