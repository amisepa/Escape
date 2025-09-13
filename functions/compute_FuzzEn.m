function FuzzEn = compute_FuzzEn(data, varargin)
% compute_FuzzEn  Computes Fuzzy Entropy (FuzzEn) across multichannel data.
%
%   FuzzEn = compute_FuzzEn(data, 'm', 2, 'n', 2, 'tau', 1, 'Kernel', 'exponential', 'Parallel', true)
%
%   Estimates signal irregularity via fuzzy similarity between embedded vectors.
%
%   Inputs:
%     - data      : EEG matrix [n_channels x n_samples]
%     - 'm'       : embedding dimension (default = 2)
%     - 'n'       : fuzziness exponent (default = 2)
%     - 'tau'     : embedding delay (default = 1)
%     - 'Kernel'  : similarity kernel ('exponential' [default] or 'gaussian')
%     - 'Parallel': true (use parfor) or false (default)
%
%   Output:
%     - FuzzEn    : Fuzzy Entropy values [n_channels x 1]
%
%   Notes:
%     - The signal is z-scored per channel internally
%     - Distance threshold r = 0.2 × std (recommended default)
%     - A fallback loop is used if vectorized computation fails (e.g., OOM)
%
%   References:
%     - Chen, W., et al., (2007). Characterization of surface EMG signal 
%       based on fuzzy entropy. IEEE Transactions on neural systems and 
%       rehabilitation engineering, 15(2).
%     - Azami, H., et al., (2019). Fuzzy entropy metrics for the analysis 
%       of biomedical signals: Assessment and comparison. IEEE Access, 7.
%
%   Copyright (C) Cedric Cannard 2025 – Escape EEGLAB Plugin (https://github.com/amisepa/Escape)

% -----------------------------
% Input parsing
% -----------------------------
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('n', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Kernel', 'exponential', @(x) ismember(lower(x), {'exponential','gaussian'}));
p.addParameter('Parallel', false, @(x) islogical(x) && isscalar(x));
p.parse(data, varargin{:});

m          = p.Results.m;
n_exp      = p.Results.n;
tau        = p.Results.tau;
kernelType = lower(p.Results.Kernel);
useParfor  = p.Results.Parallel;

% Ensure [channels x time]
if size(data,1) > size(data,2)
    data = data';
end

% -----------------------------
% Preprocessing
% -----------------------------
[nchan, ~] = size(data);
FuzzEn = nan(nchan,1);
data_z = normalize(data, 2);
r_vals = 0.2 * std(data_z, 0, 2);  % As recommended in Estevez-Baez et al. (2021)

fprintf('Computing fuzzy entropy (FuzzEn) on %g EEG channels...\n', nchan);

% -----------------------------
% Main Loop
% -----------------------------
if useParfor
    parfor iChan = 1:nchan
        FuzzEn(iChan) = compute_FuzzEn_single(data_z(iChan,:), m, r_vals(iChan), n_exp, tau, kernelType);
        fprintf('FuzzEn for channel %3d/%3d : %6.3f\n', iChan, nchan, FuzzEn(iChan));
    end
else
    progressbar('Computing FuzzEn on all channels (serial mode)')
    for iChan = 1:nchan
        FuzzEn(iChan) = compute_FuzzEn_single(data_z(iChan,:), m, r_vals(iChan), n_exp, tau, kernelType);
        progressbar(iChan / nchan)
        fprintf('FuzzEn for channel %3d/%3d : %6.3f\n', iChan, nchan, FuzzEn(iChan));
    end
end

end

% =========================================================================
function fe = compute_FuzzEn_single(signal, m, r, n, tau, kernelType)
% Helper function for per-channel fuzzy entropy
signal = signal(:);
if tau > 1, signal = signal(1:tau:end); end
N = length(signal);
p = zeros(1, 2);

for k = m:m+1
    M = N - (k - 1)*tau;
    if M <= 1
        fe = NaN;
        return
    end

    % Embedding
    X = zeros(M, k);
    for j = 1:k
        X(:, j) = signal(1 + (j - 1)*tau : N - (k - j)*tau);
    end

    try
        % Fast chebyshev distance via broadcasting
        D = max(abs(permute(X, [1 3 2]) - permute(X, [3 1 2])), [], 3);
        D(1:M+1:end) = NaN;  % Exclude diagonal

        switch kernelType
            case 'exponential'
                mu = exp(-(D.^n) / r);
            case 'gaussian'
                mu = exp(-(D.^2) / (2*r^2));
        end
        p(k - m + 1) = nanmean(mu(:));

    catch
        % Fallback: slower but memory-efficient loop
        mu = zeros(M, 1);
        for i = 1:M
            d = max(abs(X - X(i, :)), [], 2);
            d(i) = [];
            switch kernelType
                case 'exponential'
                    mu(i) = mean(exp(-(d.^n) / r));
                case 'gaussian'
                    mu(i) = mean(exp(-(d.^2) / (2*r^2)));
            end
        end
        p(k - m + 1) = mean(mu);
    end
end

if any(p == 0)
    fe = NaN;
else
    fe = log(p(1) / p(2));
end
end
