%% Simulate pink (1/f), white (0), brown (1/f²), and black noise (1/f³)
% Calculates entropy measures, and compares these signals with EEG data,
% inspired by the methodologies described in:
%   Donoghue et al. (2024). Evaluating and Comparing Measures of Aperiodic
%   Neural Activity. BioRXiv.
%
% Cedric Cannard, Jan 2025

clear; close all; clc

% Parameters
fs = 128; % Sampling frequency (Hz)
duration = 60; % Duration (seconds)
t = (0:1/fs:duration-1/fs)'; % Time vector
n = length(t); % Number of samples

% Generate Noise
white_noise = randn(n, 1); % White noise (flat spectrum)

% Pink noise
pink_noise_gen = dsp.ColoredNoise('Color', 'pink', 'SamplesPerFrame', n);
pink_noise = pink_noise_gen(); % Generate pink noise

% Brown noise
brown_noise_gen = dsp.ColoredNoise('Color', 'brown', 'SamplesPerFrame', n);
brown_noise = brown_noise_gen(); % Generate brown noise

% Generate Black Noise (1/f^3)
freqs = (0:n-1) / n; % Normalized frequency
freqs = fftshift(freqs - 0.5); % Center frequencies
freqs = ifftshift(abs(freqs) + eps); % Avoid division by zero
% spectrum = 1 ./ (freqs.*100); % Apply 1/f^3 scaling
alpha = 3;
spectrum = 1 ./ (freqs.^alpha + eps);
spectrum = spectrum .* (randn(size(spectrum)) + 1i*randn(size(spectrum))); % Add random phases
spectrum = spectrum / max(abs(spectrum)); % Normalize spectrum
black_noise = real(ifft(spectrum, 'symmetric')); % Inverse FFTblack_noise = black_noise(:); % ensure column vector
black_noise = black_noise(:); % ensure column



% Store Signals
noise_signals = [white_noise pink_noise brown_noise black_noise];
noise_labels = {'White Noise', 'Pink Noise', 'Brown Noise', 'Black Noise'};

% Time Series Plot with Vertical Offsets
figure('color','w'); hold on
offsets = [60, 30, 0, -30]; % Vertical offsets for each noise type
colors = lines(4); % Generate distinct colors for each noise type
for i = 1:4
    plot(t, noise_signals(:, i) + offsets(i), 'Color', colors(i, :), 'LineWidth', 1.5);
end
legend(noise_labels, 'Location', 'northeast');
xlabel('Time (s)');
ylabel('Amplitude (with offsets)');
title('Time Series of Different Noise Types');
grid on;
legend(noise_labels)

% Plot Power Spectra
figure('color','w')
for iSig = 1:4
    [pxx, f] = pwelch(noise_signals(:, iSig), fs, [], [], fs);
    pxx = pxx / max(pxx);
    subplot(2, 2, iSig);
    plot(f, 10*log10(pxx));
    % loglog(f, pxx)
    title(['Power Spectrum of ', noise_labels{iSig}]);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
end

% Entropy Parameters
m = 2;
r = 0.2;

% Entropy Measures using compute_ae
entropy_values = zeros(4, 1);
for iSig = 1:4
    entropy_values(iSig) = compute_ae(zscore(noise_signals(:, iSig)), m, r * std(noise_signals(:, iSig)));
end

% Display Entropy Results
disp('Entropy Values:');
for iSig = 1:4
    fprintf('%s: %f\n', noise_labels{iSig}, entropy_values(iSig));
end

% Compare with EEG Data (Placeholder: Replace with actual EEG dataset)
% Load EEG data (assumed preprocessed, e.g., filtered and downsampled)
% eeg_data = randn(length(t), 1); % Placeholder for EEG data
eeglab; close
EEG = import_muse('C:\Users\CedricCannard\Documents\eeg_data\muse1\sub-0a291a7fbc_ses-02_task-rest_run-01.csv');
for iChan = 1
    eeg_entropy(iChan,:) = compute_ae(zscore(EEG.data(iChan,:)), m, r * std(EEG.data(iChan,:)));
end

fprintf('EEG Data Entropy: %f\n', eeg_entropy);

% Plot Comparison
figure;
bar([entropy_values; eeg_entropy]);
set(gca, 'XTickLabel', [noise_labels, {'EEG Data'}]);
xlabel('Signal Type');
ylabel('Entropy Value');
title('Entropy Comparison');

