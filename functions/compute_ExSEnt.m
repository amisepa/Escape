function [HD, HA, HDA, info] = compute_ExSEnt(data, varargin)
% compute_ExSEnt  Extrema-Segmented Entropy (ExSEnt) for multichannel data.
%
%   [HD, HA, HDA, info] = compute_ExSEnt(data, 'm', 2, ...
%                                        'alpha', 0.20, 'lambda', 0.01, ...
%                                        'Parallel', false, 'Plot', false)
%
%   Inputs:
%     - data      : data matrix [n_channels x n_samples]
%     - 'm'       : embedding dimension (default = 2)
%     - 'alpha'   : scaling factor for r (r = alpha * std, default = 0.20)
%     - 'lambda'  : threshold factor for extrema detection (default = 0.01)
%     - 'Parallel': parallel computing (default = false)
%     - 'Plot'    : plot validation figure for the first channel (default = false)
%
%   Outputs:
%     - HD, HA, HDA : entropy values per channel
%     - info        : struct with per-channel diagnostics
%
%   Notes:
%     - Durations are in **samples** (not seconds).
%     - Amplitudes are net differences between extrema.
%     - Joint entropy: interleaves normalized D and A, embedding = 2*m.
%
%   Example:
%     [HD,HA,HDA,info] = compute_ExSEnt(EEG, 'm', 2, 'alpha', 0.2, ...
%                                       'lambda', 0.01, 'Plot', true);
%
%   Reference:
%     Kamali, S., Baroni, F., & Varona, P. (2025).
%     ExSEnt: Extrema-Segmented Entropy Analysis of Time Series.
%     arXiv:2509.07751
% 
%   Copyright (C) Cedric Cannard 2025 – Escape EEGLAB Plugin
%   https://github.com/amisepa/Escape

% Parse inputs
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('alpha', 0.20, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('lambda', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('Parallel', false, @(x) islogical(x) || ismember(lower(x),[false true]));
p.addParameter('Plot', false, @(x) islogical(x) || ismember(lower(x),[false true]));
p.parse(data, varargin{:});

m           = p.Results.m;
alpha       = p.Results.alpha;
lambda      = p.Results.lambda;
parallelMode= p.Results.Parallel;
doPlot      = p.Results.Plot;

if size(data,1) > size(data,2), data = data'; end
[nchan, ~] = size(data);

HD  = nan(nchan,1);
HA  = nan(nchan,1);
HDA = nan(nchan,1);
info = struct('M',nan(nchan,1),'rD',nan(nchan,1),'rA',nan(nchan,1), ...
              'rDA',nan(nchan,1),'rangeD',nan(nchan,1),'rangeA',nan(nchan,1), ...
              'segment_idx',[]);

fprintf('Computing ExSEnt on %g channels...\n', nchan);

if parallelMode
    parfor ch = 1:nchan
        [HD(ch),HA(ch),HDA(ch),info(ch)] = exsent_channel(data(ch,:), m, alpha, lambda);
    end
else
    for ch = 1:nchan
        [HD(ch),HA(ch),HDA(ch),info(ch)] = exsent_channel(data(ch,:), m, alpha, lambda);
        fprintf('Ch %3d/%3d | Segments=%4d | HD=%6.3f  HA=%6.3f  HDA=%6.3f\n', ...
                ch, nchan, info(ch).M, HD(ch), HA(ch), HDA(ch));
    end
end

% ===================== Validation plot (first channel) ==================== %
if doPlot
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
    yyaxis left; plot(D_vals,'b-o','LineWidth',1.2); ylabel('Durations (samples)');
    yyaxis right; plot(A_vals,'r-s','LineWidth',1.2); ylabel('Amplitudes (Δ)');
    xlabel('Segment index');
    title('Duration and amplitude sequences (D and A)');
    grid on;
end

end



% ========================= Channel-level helper ========================= %
function [HD,HA,HDA,info] = exsent_channel(x, m, alpha, lambda)
% Segment signal into durations/amplitudes
[D_vals, A_vals, segment_idx] = extract_DA(x, lambda);
M = numel(D_vals);

if M <= m+1
    [HD,HA,HDA] = deal(NaN);
    info = struct('M',M,'rD',NaN,'rA',NaN,'rDA',NaN, ...
                  'rangeD',NaN,'rangeA',NaN,'segment_idx',segment_idx);
    return;
end

% Tolerances
rD = alpha * std(D_vals);
rA = alpha * std(A_vals);

% Use your compute_SampEn_single backend
HD = compute_SampEn_single(D_vals, m, rD, 1);
HA = compute_SampEn_single(A_vals, m, rA, 1);

% Normalize then interleave [D_norm, A_norm]
D_norm = normalize(D_vals);
A_norm = normalize(A_vals);
joint_data = reshape([D_norm(:) A_norm(:)]', [], 1);

% Joint entropy with embedding = 2*m
HDA = compute_SampEn_single(joint_data, 2*m, alpha, 1);

% Diagnostics
info = struct('M',M,'rD',rD,'rA',rA,'rDA',alpha, ...
              'rangeD',range(D_vals),'rangeA',range(A_vals), ...
              'segment_idx',segment_idx);
end


function [D_vals, A_vals, segment_idx] = extract_DA(signal, lambda)
% extract_DA  Extract segment durations and amplitudes based on extrema.
%
%   [D_vals, A_vals, segment_idx] = extract_DA(signal, lambda)
%
%   Inputs:
%     - signal : 1D time series vector (row or column)
%     - lambda : scaling factor for thresholding small derivatives
%                (default = 0.001). Threshold = lambda * IQR(diff(signal)).
%
%   Outputs:
%     - D_vals      : array of segment durations (in samples)
%     - A_vals      : array of segment net amplitude changes
%     - segment_idx : indices of detected extrema (endpoints of segments)
%
%   Notes:
%     - A segment is defined as a monotone run of samples between two extrema.
%     - Extrema are detected at sign changes of the derivative that exceed
%       the noise threshold.
%     - Durations are in samples. To convert to time, divide by fs outside.
%
%   Example:
%     [D,A,idx] = extract_DA(x, 0.01);
%     plot(x); hold on; plot(idx, x(idx),'ro');
%
%   Reference: Kamali, S., et al., (2025). ExSEnt: Extrema-Segmented Entropy
%              Analysis of Time Series. arXiv:2509.07751
% 
% Copyright (C) Cedric Cannard 2025 – Escape EEGLAB Plugin (https://github.com/amisepa/Escape)

if nargin < 2
    lambda = 0.001;
end

% Ensure row vector
signal = signal(:)';

% Step 1: Raw derivative
dx = diff(signal);

% Step 2: Adaptive threshold from IQR
threshold = lambda * iqr(dx);

% Step 3: Segment at sign changes with sufficient magnitude
D_vals = [];
A_vals = [];
segment_idx = [];

start_idx = 1;
n = length(signal);

for i = 2:length(dx)
    if sign(dx(i)) ~= sign(dx(i-1)) && abs(dx(i) - dx(i-1)) > threshold
        % End current segment at i
        duration = i - start_idx;
        amplitude = signal(i) - signal(start_idx);

        D_vals(end+1) = duration; %#ok<*AGROW>
        A_vals(end+1) = amplitude;
        segment_idx(end+1) = i;

        % Start new segment
        start_idx = i;
    end
end
end

%% helpers

function SampEn = compute_SampEn_single(signal, m, r, tau)
signal = signal(:)';

try
    N = length(signal);
    if N <= m + 1, SampEn = NaN; return; end

    X_m  = buffer(signal, m,  m-1, 'nodelay')';
    X_m1 = buffer(signal, m+1, m,   'nodelay')';

    D_m  = pdist(X_m,  'chebychev');
    D_m1 = pdist(X_m1, 'chebychev');

    A = mean(D_m  < r);
    B = mean(D_m1 < r);

    if A > 0 && B > 0
        SampEn = -log(B / A);
    else
        SampEn = NaN;
    end

catch
    warning('Optimized SampEn failed. Trying fallback method.');
    SampEn = fallback_SampEn(signal, m, r, tau);
end
end

function SampEn = fallback_SampEn(signal, m, r, tau)
if tau > 1, signal = signal(1:tau:end); end
n = length(signal);
p = zeros(1,2);
sMat = zeros(m+1,n-m);

for i = 1:m+1
    sMat(i,:) = signal(i:n-m+i-1);
end

for k = m:m+1
    count = zeros(1,n-m);
    tempMat = sMat(1:k,:);
    for i = 1:n-k
        dist = max(abs(tempMat(:,i+1:n-m) - repmat(tempMat(:,i),1,n-m-i)));
        count(i) = sum(dist < r)/(n-m);
    end
    p(k-m+1) = sum(count)/(n-m);
end

if any(p == 0)
    SampEn = NaN;
else
    SampEn = log(p(1)/p(2));
end
end