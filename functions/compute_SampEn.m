function SampEn = compute_SampEn(data, varargin)
% compute_SampEn  Computes Sample Entropy (SampEn) across multichannel data.
%
%   SampEn = compute_SampEn(data, 'm', 2, 'tau', 1, 'Parallel', false)
%
%   Inputs:
%     - data      : EEG data matrix [n_channels x n_samples]
%     - 'm'       : embedding dimension (default = 2)
%     - 'tau'     : time lag (default = 1; used in fallback only)
%     - 'Parallel': parallel computing turned ON (true; default) or not (false)
%
%   Output:
%     - SampEn : Sample Entropy vector [n_channels x 1]
%
%   Improvements:
%     - Supports multichannel input (vectorized where possible)
%     - Computes r and z-score internally per channel
%     - Automatically handles parfor for large datasets
%     - Tracks progress using command line + progressbar
%
%   Copyright (C) Cedric Cannard, 2022–2025
%   Part of the ASCENT EEGLAB plugin (https://github.com/cedriccannard/ascent)


% Parse inputs
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ndims(x) == 2);
p.addParameter('m', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tau', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Parallel', false, @(x) ismember(lower(x), [false,true]));
p.parse(data, varargin{:});

m = p.Results.m;
tau = p.Results.tau;
parallelMode = p.Results.Parallel;

if size(data,1) > size(data,2)
    data = data';  % [n_channels x n_samples]
end

[nchan, ~] = size(data);
SampEn = nan(nchan, 1);
data_z = normalize(data, 2);             % z-score across time
r_vals = 0.15 * std(data_z, 0, 2);        % per channel

fprintf('Computing sample entropy (SampEn) on %g EEG channels...\n', nchan);

if parallelMode
    parfor iChan = 1:nchan
        SampEn(iChan) = compute_SampEn_single(data_z(iChan,:), m, r_vals(iChan), tau);
        fprintf('SampEn for channel %3d/%3d : %6.3f\n', iChan, nchan, SampEn(iChan));  
    end
else
    progressbar('Computing SampEn on all channels (serial mode)')
    for iChan = 1:nchan
        SampEn(iChan) = compute_SampEn_single(data_z(iChan,:), m, r_vals(iChan), tau);
        progressbar(iChan / nchan)
        fprintf('SampEn for channel %3d/%3d : %6.3f\n', iChan, nchan, SampEn(iChan));
    end
end
end

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
