function [HD, HA, H_joint,M,range_D,range_A,r_D,r_A, segment_ids] = compute_ExSEnt(signal,lambda,m,r)

% COMPUTE_ExSEnt
%   Computes Sample Entropy (SampEn) on:
%     (1) durations between successive extrema segments, D
%     (2) cumulative amplitudes across those segments, A
%     (3) the joint (D,A) sequence
%
% INPUTS
%   signal   : 1D time series (row or column vector)
%   lambda   : threshold parameter used inside extract_DA (as % of IQR) to
%              enforce noise tolerance in extrema segmentation
%   m        : embedding dimension for SampEn on the univariate series
%   r        : scaling factor for the tolerance r (i.e., r = r * std)
%
% OUTPUTS
%   HD         : SampEn of durations D at embedding m (scalar)
%   HA         : SampEn of cumulative amplitudes A at embedding m (scalar)
%   H_joint    : SampEn of the concatenated, normalized joint series (D,A)
%                at embedding 2*m (scalar)
%   M          : number of detected extrema segments
%   range_D    : range(D) = max(D) - min(D)
%   range_A    : range(A) = max(A) - min(A)
%   r_D, r_A   : SampEn tolerances used for D and A respectively
%   segment_ids: indices/labels of segments returned by extract_DA
%
% DEPENDENCIES (must exist on path)
%   extract_DA.m      → [D_vals, A_vals, segment_ids] = extract_DA(signal, lambda)
%   sample_entropy.m  → H = sample_entropy(x, m, r)
%
% NOTES
%   • r choice: standard practice is r = r * std(x).
%   • Joint entropy: we embed a single 1D sequence built by interleaving
%     normalized D and A, then use mAD = 2*m to match total state size.
%   • Normalization: To make r_joint interpretable across datasets, we put
%     D and A on comparable scales before interleaving.
%     IMPORTANT: MATLAB's normalize(x) by default does z-score (mean 0, std 1).
%     If you intend [0,1] range normalization, use normalize(x,'range').

% Extract durations (D_vals) and cumulative amplitudes (A_vals)
[D_vals, A_vals, segment_ids] = extract_DA(signal, lambda);

% Check if valid segments were found
if isempty(D_vals)
    error('No valid segments found in the raw signal with the current threshold.');
end

% Set parameters for sample entropy computation
r_D = r * std(D_vals); % tolerance for durations (D)
r_A = r * std(A_vals); % tolerance for amplitudes (A)

% Compute sample entropy for durations (D)
i=1;
% for m=2:20
HD(i) =  sample_entropy(D_vals, m, r_D);%sample_entropy(D_vals, m, r_D);

% Compute sample entropy for cumulative amplitudes (A)
HA(i) = sample_entropy(A_vals, m, r_A);%sampen


% Compute oint sample entropy for the normalized paired (D, A)
D_vals_norm=normalize(D_vals);  % Normalizes to [0,1]
A_vals_norm=normalize(A_vals);  % Normalizes to [0,1]
joint_data_norm = reshape([D_vals_norm(:) A_vals_norm(:)].',[],1);
mAD = 2* m;

% r_joint = r * std(joint_data(:)); % overall tolerance for joint data
r_joint = r; % Tolerance for joint data since (std=1 for normalized data)
H_joint(i) = sample_entropy(joint_data_norm, mAD, r_joint);
% i=i+1;end
% set(gca,'fontsize',12)
% figure;subplot(311);plot(2:20,HD,'r-s','LineWidth',2);ylabel('H_D');hold on;subplot(312);plot(2:20,HA,'b-v','LineWidth',2);ylabel('H_A');
% subplot(313);plot(2:20,H_joint,'k-d','LineWidth',2);ylabel('H_{DA}');xlabel('m');

M=length(D_vals); %number of segments
range_D = range(D_vals);%range of duration values
range_A = range(A_vals);%range of amplitudes


end


% extract_DA extracts duration (D) and cumulative amplitude (A) pairs from a
% 1D time series by segmenting based on sign changes in the derivative.

function [duration_sequence, amp_sequence,segments_idx] = extract_DA(time_series, lambda)

% INPUT:
%   time_series - 1D array representing the time series data.
%   lambda      - Scaling factor for thresholding.
%
% OUTPUT:
%   duration_sequence   - Array containing the duration (number of points) of each segment.
%   amp_sequence     - Array containing the cumulative amplitude (sum of values) of each segment.
if nargin<2
    lambda=0.001;
end
% Step 1: Compute the raw derivative (no smoothing)
derivative = diff(time_series);

% Step 2: Compute an adaptive threshold using robust quantiles
% This threshold is based on the difference between the 75th and 25th percentiles of diff(x)
threshold = lambda * iqr(derivative);

% Step 3: Segment the signal based on sign changes or low derivative values.
% No minimum segment length is enforced.
duration_sequence = [];
amp_sequence = [];
segments_idx =[];
start_idx = 1;

for i = 2:length(derivative)
    if sign(derivative(i)) ~= sign(derivative(i -1)) && (abs(derivative(i) - derivative(i-1))>threshold)
        % End the current segment
        duration = i - start_idx;
        amplitude = time_series(i) - time_series(start_idx);
        duration_sequence = [duration_sequence, duration];
        amp_sequence = [amp_sequence, amplitude];
        segments_idx = [segments_idx, i ];
        start_idx = i;  % Start a new segment
    end
end
end

% Computes sample entropy for a 1D data series.

function samp_ent = sample_entropy(data, m, r)

N = length(data);
if N <= m+1
    samp_ent = NaN;
    return;
end

% Construct templates of length m.
X = zeros(N-m+1, m);
for i = 1:(N-m+1)
    X(i,:) = data(i:i+m-1);
end

% Count similar template pairs for m
B = 0;
for i = 1:size(X,1)
    for j = i+1:size(X,1)
        if max(abs(X(i,:) - X(j,:))) <= r
            B = B + 1;
        end
    end
end

% Construct templates of length m+1.
X1 = zeros(N-m, m+1);
for i = 1:(N-m)
    X1(i,:) = data(i:i+m);
end

% Count similar template pairs for m+1.
A = 0;
for i = 1:size(X1,1)
    for j = i+1:size(X1,1)
        if max(abs(X1(i,:) - X1(j,:))) <= r
            A = A + 1;
        end
    end
end

if B == 0 || A == 0
    samp_ent = Inf;
else
    samp_ent = -log(A / B);
end

end