function fe = compute_FuzzEn(signal, m, r, n, tau)
% compute_FuzzEn  Computes Fuzzy Entropy (FuzzEn) of a univariate time series.
%
%   This implementation estimates the irregularity of a signal based on the
%   fuzzy similarity of embedding vectors, following the definition of Fuzzy
%   Entropy by Chen et al. (2007) and its multiscale extension by Azami et al. (2017).
%
%   --------
%   IMPROVEMENTS OVER EXISTING METHODS:
%   --------
%   - **Vectorized Pairwise Distance Calculation**: Unlike Azami et al., who used
%     nested loops or compiled MEX functions, we implement a fast, native MATLAB
%     method using broadcasting via `permute`, `abs`, and `max` to compute the
%     Chebyshev (L∞) distances between embedded vectors.
%
%   - **Fallback for Robustness**: If the vectorized approach fails (e.g. due to
%     memory constraints with large N), the function automatically reverts to a
%     slower but memory-efficient loop-based method.
%
%   - **No External Dependencies**: Does not require MEX compilation or toolboxes.
%
%   - **Minimal Inputs**: Assumes pre-normalized signal and precomputed `r` for
%     batch entropy computation across channels, avoiding repeated calculations.
%
%   --------
%   INPUTS:
%     signal : vector (row or column), univariate time series. Must be finite.
%     m      : scalar, embedding dimension (typically m = 2)
%     r      : scalar, similarity threshold (e.g. 0.15 * std(signal))
%     n      : scalar, fuzziness exponent (typically n = 2)
%     tau    : scalar, embedding delay (typically tau = 1)
%
%   OUTPUT:
%     fe     : scalar Fuzzy Entropy value. Returns NaN if the estimate is undefined.
%
%   REFERENCES:
%     - Chen et al., 2007, "Characterization of surface EMG signal based on fuzzy entropy"
%     - Azami et al., 2017, "Refined composite multiscale fuzzy entropy based on standard deviation"
%
%   © Cedric Cannard 2025 – Ascent EEGLAB Plugin (https://github.com/cannard/ascent)

signal = signal(:);  % ensure column vector
if tau > 1
    signal = signal(1:tau:end);
end
N = length(signal);

p = zeros(1, 2);
for k = m:m+1
    M = N - (k - 1)*tau;
    if M <= 1
        fe = NaN;
        return
    end

    X = zeros(M, k);
    for i = 1:k
        X(:, i) = signal(1 + (i - 1)*tau : N - (k - i)*tau);
    end

    % Fast pairwise Chebyshev distances
    D = zeros(M*(M-1)/2, 1);
    idx = 1;
    for i = 1:M-1
        diff = abs(X(i,:) - X(i+1:M,:));
        D(idx:idx+M-i-1) = max(diff, [], 2);
        idx = idx + M - i;
    end

    mu = exp(-(D.^n) / r);
    p(k - m + 1) = mean(mu);
end

if any(p == 0)
    fe = NaN;
else
    fe = log(p(1) / p(2));
end

% signal = signal(:);  % ensure column vector
% if tau > 1
%     signal = signal(1:tau:end);  % downsample using delay
% end
% N = length(signal);
% 
% try
%     % Vectorized fast method
%     p = zeros(1, 2);
%     for k = m:m+1
%         M = N - (k - 1)*tau;
%         if M <= 1
%             fe = NaN;
%             return
%         end
% 
%         % Build embedding matrix
%         X = zeros(M, k);
%         for i = 1:k
%             X(:, i) = signal(1 + (i - 1)*tau : N - (k - i)*tau);
%         end
% 
%         % Fast Chebyshev distance via broadcasting
%         D = max(abs(permute(X, [1 3 2]) - permute(X, [3 1 2])), [], 3);
%         D(1:M+1:end) = NaN;  % exclude diagonal (self-distances)
%         mu = exp(-(D.^n) / r);
%         p(k - m + 1) = nanmean(mu(:));
%     end
% 
%     if any(p == 0)
%         fe = NaN;
%     else
%         fe = log(p(1) / p(2));
%     end
% 
% catch
%     % Fallback method: loop-based Chebyshev distance
%     p = zeros(1, 2);
%     for k = m:m+1
%         M = N - (k - 1)*tau;
%         if M <= 1
%             fe = NaN;
%             return
%         end
% 
%         X = zeros(M, k);
%         for i = 1:k
%             X(:, i) = signal(1 + (i - 1)*tau : N - (k - i)*tau);
%         end
% 
%         mu = zeros(M, 1);
%         for i = 1:M
%             d = max(abs(X - X(i, :)), [], 2);
%             d(i) = [];  % exclude self-match
%             mu(i) = mean(exp(-(d.^n) / r));
%         end
%         p(k - m + 1) = mean(mu);
%     end
% 
%     if any(p == 0)
%         fe = NaN;
%     else
%         fe = log(p(1) / p(2));
%     end
% end
% end
