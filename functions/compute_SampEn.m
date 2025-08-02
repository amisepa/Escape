%% compute_SampEn.m
% Computes sample entropy (SampEn) of a univariate signal.
%
% Sample entropy quantifies the negative logarithm of the conditional probability 
% that two sequences of length m that match within a tolerance r 
% will also match at the next point (length m+1), excluding self-matches.
%
% Inputs:
%   signal - univariate signal; a vector of size 1 x N or N x 1
%   m      - embedding dimension (default = 2)
%   r      - tolerance (commonly 0.1–0.25 * std(signal))
%   tau    - time lag (optional; used only in fallback method; default = 1)
%
% Outputs:
%   se     - sample entropy
%
% Dependencies:
%   NEW method requires:
%     - pdist (Statistics & Machine Learning Toolbox)
%     - buffer (Signal Processing Toolbox)
%
% If unavailable, it automatically falls back to older methods.
%
% Copyright (C) Cedric Cannard, August 2022 – July 2025

function se = compute_SampEn(signal, m, r, tau)

% Default parameters
if nargin < 2, m = 2; end
if nargin < 3, r = 0.15 * std(signal); end
if nargin < 4, tau = 1; end

% Ensure row vector
if size(signal, 1) > size(signal, 2)
    signal = signal';
end

% Try NEW optimized method (fastest)
try
    % Check signal length
    N = length(signal);
    if N <= m + 1
        error('Signal too short for given embedding dimension.');
    end

    % Construct embedded vectors
    X_m  = buffer(signal, m,  m-1, 'nodelay')';
    X_m1 = buffer(signal, m+1, m,   'nodelay')';

    % Pairwise Chebyshev distances
    D_m  = pdist(X_m,  'chebychev');
    D_m1 = pdist(X_m1, 'chebychev');

    % Match proportions
    A = mean(D_m  < r);
    B = mean(D_m1 < r);

    % Entropy
    if A == 0 || B == 0
        se = NaN;
    else
        se = -log(B / A);
    end

catch
    warning('Fast SampEn method failed. Trying alternative method (Vest et al.)...');
    try
        % Fallback to VEST ET AL. fast method
        if size(signal,1) > size(signal,2), signal = signal'; end
        xx = convert_to_lagged_form(signal, m)';
        Dxx = pdist(xx,'chebychev');
        yy = convert_to_lagged_form(signal, m+1)';
        Dyy = pdist(yy,'chebychev');
        A = mean(Dxx < r);
        B = mean(Dyy < r);
        se = -log(B / A);
    catch
        warning('Vest method also failed. Reverting to original method...');
        % Original method with no toolbox dependencies
        if tau > 1, signal = downsamp(signal, tau); end
        n = length(signal);
        p = zeros(1,2);
        sMat = zeros(m+1,n-m);
        parfor i = 1:m+1
            sMat(i,:) = signal(i:n-m+i-1);
        end
        for k = m:m+1
            count = zeros(1,n-m);
            tempMat = sMat(1:k,:);
            parfor i = 1:n-k
                dist = max(abs(tempMat(:,i+1:n-m) - repmat(tempMat(:,i),1,n-m-i)));
                count(i) = sum(dist < r)/(n-m);
            end
            p(k-m+1) = sum(count)/(n-m);
        end
        se = log(p(1)/p(2));
    end
end
end

%% Helper functions (used in fallback methods)
function yy = convert_to_lagged_form(y, k)
[s, T] = size(y);
bs = s*ones(1,k);
yy = zeros(k*s, T-k+1);
for i = 1:k
    yy(block(i,bs), :) = y(:, k-i+1:end-i+1); 
end
end

function sub = block(blocks, block_sizes)
blocks = blocks(:)';
block_sizes = block_sizes(:)';
skip = [0 cumsum(block_sizes)];
start = skip(blocks)+1;
fin = start + block_sizes(blocks) - 1;
sub = [];
for j = 1:length(blocks)
    sub = [sub start(j):fin(j)];
end
end

function y = downsamp(x, n, phase)
if nargin < 3, phase = 0; end
if isvector(x)
    y = x(phase + 1:n:end);
else
    y = x(phase + 1:n:end,:);
end
end
