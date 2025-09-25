clear; close all; clc

%% STREAMLINED ENTROPY ALGORITHM VALIDATION SCRIPT
% Validates new MSE algorithm against theoretical expectations from synthetic signals
% No Costa comparison - focuses on intrinsic algorithm behavior

%% Configuration
config = struct();
config.fs = 128;                    % Hz
config.dur = 180;                   % seconds (3 min) - main signals
config.tolerance_threshold = 0.10;  % for flagging large deviations
config.alpha = 0.05;               % significance level for statistical tests
config.n_bootstrap = 1000;         % bootstrap samples for confidence intervals

% Algorithm parameters to test
config.r_values = [0.10, 0.15, 0.20, 0.25];  % r parameters to test
config.m_values = [2, 3];                     % embedding dimensions to test
config.default_r = 0.15;                      % default r for main analysis
config.default_m = 2;                         % default m for main analysis
config.default_nscales = 30;
config.default_coarsing = 'mean';    % 'mean' or 'std'
config.default_filtmode = 'none';    % 'none' or 'narrowband' (Kosciessa style)

% Edge case testing
config.short_dur = [10, 30, 60];              % short signal durations (seconds)

fprintf('=== Streamlined Entropy Algorithm Validation ===\n\n');

%% Generate Signal Suite with Known Properties
fprintf('Generating signal suite with known entropy properties...\n');
N = config.fs * config.dur;
rng(1) % reproducibility

sig = struct();

% 1) Pure periodic (very regular, lowest entropy)
sig.periodic = sin(2*pi*10*(0:N-1)/config.fs);

% 2) White Gaussian noise (moderate entropy, decreasing with scale)
sig.white = randn(1,N);

% 3) Pink (1/f) noise (persistent, moderate entropy)
K = floor(N/2);
mag = (1:K).^(-1/2);        % amplitude ~ f^{-beta/2}, beta=1
mag(1) = 1;                 % guard DC
phi = 2*pi*rand(1,K);
pos = mag .* exp(1j*phi);
X = zeros(1,N);
X(1) = 0;                   % remove DC
X(2:K+1) = pos;
if mod(N,2)==0
    X(K+1) = real(X(K+1));       % Nyquist real
    X(K+2:end) = conj(pos(K-1:-1:1));
else
    X(K+2:end) = conj(pos(K:-1:1));
end
sig.pink = zscore(real(ifft(X))*sqrt(N));

% 4) Brown noise (integrated white; increasing entropy with scale)
sig.brown = zscore(cumsum(randn(1,N)));

% 5) AR(1) processes: weak (phi=0.3) and strong (phi=0.9) correlation
phi = 0.3;
e = randn(1,N);
x = zeros(1,N);
x(1) = randn * (1/sqrt(1-phi^2));
for t = 2:N, x(t) = phi*x(t-1) + e(t); end
sig.AR03 = zscore(x);

phi = 0.9;
e = randn(1,N);
x = zeros(1,N);
x(1) = randn * (1/sqrt(1-phi^2));
for t = 2:N, x(t) = phi*x(t-1) + e(t); end
sig.AR09 = zscore(x);

% 6) Sinusoid + small noise (quasi-periodic)
sig.sinNoise = sin(2*pi*8*(0:N-1)/config.fs) + 0.2*randn(1,N);

% 7) Logistic-map chaos (complex, high entropy)
r = 4.0;
burn = 200;
x0 = 0.5 + 1e-3*randn;            % non-degenerate
x0 = min(max(x0,1e-6), 1-1e-6);
x = zeros(1,N+burn); x(1) = x0;
for t = 1:N+burn-1, x(t+1) = r*x(t)*(1 - x(t)); end
sig.logistic = zscore(x(burn+1:end));

% 8) Blue noise (anti-persistent; beta = -1)
K = floor(N/2);
mag = (1:K).^(+1/2);        % amplitude ~ f^{-beta/2}, beta=-1 -> f^{+1/2}
mag(1) = 1;
phi = 2*pi*rand(1,K);
pos = mag .* exp(1j*phi);
X = zeros(1,N);
X(1) = 0;
X(2:K+1) = pos;
if mod(N,2)==0
    X(K+1) = real(X(K+1));
    X(K+2:end) = conj(pos(K-1:-1:1));
else
    X(K+2:end) = conj(pos(K:-1:1));
end
sig.blue = zscore(real(ifft(X))*sqrt(N));

names = fieldnames(sig);
num_sigs = numel(names);

% Expected theoretical properties
expected_ranking = {'periodic', 'brown', 'AR09', 'pink', 'AR03', 'sinNoise', 'blue', 'logistic', 'white'};
scale_behaviors = containers.Map(...
    {'periodic', 'white', 'pink', 'brown', 'AR03', 'AR09', 'sinNoise', 'logistic', 'blue'}, ...
    {'flat_low', 'decreasing', 'mild_decrease', 'increasing', 'decreasing', 'decreasing', 'hump', 'decreasing', 'decreasing'});

fprintf('Expected SampEn ranking (low→high): %s\n\n', strjoin(expected_ranking, ' < '));

%% Signal Characterization and Visualization
fprintf('Analyzing signal properties...\n');

% Create comprehensive signal overview
figure('Color','w','Name','Signal Suite Analysis','Position',[50 50 1400 900]);

% Signal time series plots (3x3 grid)
for i = 1:num_sigs
    subplot(3,3,i);
    x = zscore(sig.(names{i}));
    plot(x(1:min(1000,end))); % Show first 1000 samples for clarity
    title(names{i},'Interpreter','none');
    grid on; axis tight;
end

% Create second figure for spectral analysis
figure('Color','w','Name','Signal Spectral Analysis','Position',[100 100 1200 600]);

% Power spectral densities
subplot(1,2,1);
hold on;
colors = lines(num_sigs);
spectral_stats = struct();
for i = 1:num_sigs
    x = sig.(names{i}) - mean(sig.(names{i}));
    [Pxx,f] = pwelch(x, 2048, 1024, 2048, config.fs);
    loglog(f(2:end), Pxx(2:end), 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Store spectral slope for analysis
    fband = f(f >= 1 & f <= 30);
    Pband = Pxx(f >= 1 & f <= 30);
    if length(fband) > 10
        p = polyfit(log10(fband), log10(Pband), 1);
        spectral_stats.(names{i}).slope = p(1);
    else
        spectral_stats.(names{i}).slope = NaN;
    end
end
grid on; xlabel('Frequency (Hz)'); ylabel('PSD');
title('Power Spectral Densities'); 
legend(names, 'Interpreter','none','Location','northeast');

% Autocorrelation functions
subplot(1,2,2);
hold on;
max_lag = min(100, floor(N/20));
for i = 1:num_sigs
    x = sig.(names{i}) - mean(sig.(names{i}));
    [acf,lags] = xcorr(x, max_lag, 'normalized');
    plot(lags, acf, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Store lag-1 autocorrelation
    lag1_idx = find(lags == 1);
    if ~isempty(lag1_idx)
        spectral_stats.(names{i}).lag1 = acf(lag1_idx);
    else
        spectral_stats.(names{i}).lag1 = NaN;
    end
end
grid on; xlabel('Lag'); ylabel('Autocorrelation');
title('Autocorrelation Functions');
legend(names, 'Interpreter','none','Location','northeast');

% Print spectral analysis results
fprintf('\nSpectral Analysis Results:\n');
fprintf('%-12s %8s %8s %12s\n', 'Signal', 'Slope', 'Lag1', 'Classification');
fprintf(repmat('-', 1, 50)); fprintf('\n');
for i = 1:num_sigs
    slope = spectral_stats.(names{i}).slope;
    lag1 = spectral_stats.(names{i}).lag1;
    
    % Classify signal type
    if abs(lag1) < 0.1 && abs(slope) < 0.5
        class = 'white-like';
    elseif lag1 > 0.8
        class = 'persistent';
    elseif lag1 < -0.3
        class = 'antipersist';
    elseif slope < -1.5
        class = '1/f-like';
    elseif abs(slope) < 0.1
        class = 'broadband';
    else
        class = 'other';
    end
    
    if isnan(slope)
        slope_str = '    N/A';
    else
        slope_str = sprintf('%+7.2f', slope);
    end
    if isnan(lag1)
        lag1_str = '   N/A';
    else
        lag1_str = sprintf('%+7.2f', lag1);
    end
    fprintf('%-12s %8s %8s %12s\n', names{i}, slope_str, lag1_str, class);
end
fprintf('\n');

%% Generate Edge Cases
fprintf('Generating edge case signals...\n');
edge_signals = struct();

% Short signals for robustness testing
for i = 1:length(config.short_dur)
    dur = config.short_dur(i);
    N_short = config.fs * dur;
    edge_signals.(sprintf('short_%ds', dur)) = struct();
    edge_signals.(sprintf('short_%ds', dur)).white = randn(1,N_short);
    edge_signals.(sprintf('short_%ds', dur)).periodic = sin(2*pi*10*(0:N_short-1)/config.fs);
end

% Different amplitude distributions
N_edge = config.fs * 60; % 1 minute
edge_signals.uniform = 2*(rand(1,N_edge)-0.5)*sqrt(3); % same variance as N(0,1)
edge_signals.exponential = zscore(exprnd(1, 1, N_edge));
edge_signals.bimodal = zscore([randn(1,floor(N_edge/2))-2, randn(1,ceil(N_edge/2))+2]);

%% Main MSE Computation
fprintf('Computing MSE for all signals...\n');

mse_params = struct('m', config.default_m, 'r', config.default_r, 'tau', 1, ...
                   'coarsing', config.default_coarsing, ...
                   'num_scales', config.default_nscales, ...
                   'filter_mode', config.default_filtmode);

% Main signal analysis
num_scales = mse_params.num_scales;
MSE_results = nan(num_sigs, num_scales);
computation_times = zeros(num_sigs, 1);

fprintf('Processing main signals:\n');
for i = 1:num_sigs
    fprintf('  %s... ', names{i});
    tic;
    try
        MSE_results(i, :) = compute_mse(sig.(names{i}), 'Fs', config.fs, 'm', mse_params.m, ...
            'tau', mse_params.tau, 'r', mse_params.r, 'num_scales', num_scales, ...
            'coarsing', mse_params.coarsing, 'filter_mode', mse_params.filter_mode, ...
            'Parallel', false, 'Progress', false);
        computation_times(i) = toc;
        fprintf('%.2fs\n', computation_times(i));
    catch ME
        computation_times(i) = toc;
        fprintf('FAILED: %s\n', ME.message);
        MSE_results(i, :) = NaN;
    end
end

%% Parameter Sensitivity Analysis
fprintf('\nTesting parameter sensitivity...\n');

% Test r parameter sensitivity on subset of signals
test_signals = {'white', 'periodic', 'brown'};
r_sensitivity = struct();

fprintf('r parameter sensitivity:\n');
for r_idx = 1:length(config.r_values)
    r_val = config.r_values(r_idx);
    fprintf('  r = %.2f: ', r_val);
    r_key = sprintf('r_%03d', round(r_val*100));
    r_sensitivity.(r_key) = struct();
    
    for s_idx = 1:length(test_signals)
        sig_name = test_signals{s_idx};
        sig_idx = find(strcmpi(names, sig_name));
        if ~isempty(sig_idx)
            try
                mse_vals = compute_mse(sig.(sig_name), 'Fs', config.fs, 'm', mse_params.m, ...
                    'r', r_val, 'num_scales', config.default_nscales, 'coarsing', mse_params.coarsing, ...
                    'filter_mode', mse_params.filter_mode, 'Parallel', false, 'Progress', false);
                r_sensitivity.(r_key).(sig_name) = mse_vals;
            catch
                r_sensitivity.(r_key).(sig_name) = NaN(1,config.default_nscales);
            end
        end
    end
    fprintf('OK\n');
end

% Test m parameter sensitivity
m_sensitivity = struct();
fprintf('m parameter sensitivity:\n');
for m_idx = 1:length(config.m_values)
    m_val = config.m_values(m_idx);
    fprintf('  m = %d: ', m_val);
    m_key = sprintf('m_%d', m_val);
    m_sensitivity.(m_key) = struct();
    
    for s_idx = 1:length(test_signals)
        sig_name = test_signals{s_idx};
        if isfield(sig, sig_name)
            try
                mse_vals = compute_mse(sig.(sig_name), 'Fs', config.fs, 'm', m_val, ...
                    'r', mse_params.r, 'num_scales', config.default_nscales, 'coarsing', mse_params.coarsing, ...
                    'filter_mode', mse_params.filter_mode, 'Parallel', false, 'Progress', false);
                m_sensitivity.(m_key).(sig_name) = mse_vals;
            catch
                m_sensitivity.(m_key).(sig_name) = NaN(1,config.default_nscales);
            end
        end
    end
    fprintf('OK\n');
end

%% Edge Case Testing
fprintf('\nTesting edge cases:\n');
edge_results = struct();
edge_names = fieldnames(edge_signals);

for i = 1:numel(edge_names)
    edge_name = edge_names{i};
    fprintf('  %s: ', edge_name);
    
    if isstruct(edge_signals.(edge_name))
        % Multiple signals in this edge case
        sub_names = fieldnames(edge_signals.(edge_name));
        edge_results.(edge_name) = struct();
        success_count = 0;
        for j = 1:numel(sub_names)
            x = edge_signals.(edge_name).(sub_names{j});
            try
                num_scales_edge = min(num_scales, floor(length(x)/10));
                mse_vals = compute_mse(x, 'Fs', config.fs, 'm', mse_params.m, ...
                    'r', mse_params.r, 'num_scales', num_scales_edge, ...
                    'coarsing', mse_params.coarsing, 'filter_mode', mse_params.filter_mode, ...
                    'Parallel', false, 'Progress', false);
                edge_results.(edge_name).(sub_names{j}) = mse_vals;
                success_count = success_count + 1;
            catch
                edge_results.(edge_name).(sub_names{j}) = NaN;
            end
        end
        fprintf('%d/%d OK\n', success_count, numel(sub_names));
    else
        % Single signal
        x = edge_signals.(edge_name);
        try
            num_scales_edge = min(num_scales, floor(length(x)/10));
            mse_vals = compute_mse(x, 'Fs', config.fs, 'm', mse_params.m, ...
                'r', mse_params.r, 'num_scales', num_scales_edge, ...
                'coarsing', mse_params.coarsing, 'filter_mode', mse_params.filter_mode, ...
                'Parallel', false, 'Progress', false);
            edge_results.(edge_name) = mse_vals;
            fprintf('OK\n');
        catch ME
            edge_results.(edge_name) = NaN;
            fprintf('FAILED\n');
        end
    end
end

%% Comprehensive Results Visualization
fprintf('\nCreating comprehensive visualizations...\n');

% Main results plot
figure('Color','w','Name','MSE Algorithm Validation Results','Position',[100 100 1600 1000]);

% MSE curves for all signals
subplot(2,4,1);
scales = 1:num_scales;
hold on;
for i = 1:num_sigs
    plot(scales, MSE_results(i,:), 'o-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', names{i});
end
xlabel('Scale'); ylabel('SampEn'); title('MSE Curves');
legend('Location', 'eastoutside', 'Interpreter', 'none'); grid on;

% Scale 1 entropy values (ranked)
subplot(2,4,2);
scale1_values = MSE_results(:, 1);
valid_idx = ~isnan(scale1_values);
valid_vals = scale1_values(valid_idx);
[sorted_vals, temp_idx] = sort(valid_vals);
sorted_idx = find(valid_idx);
sorted_idx = sorted_idx(temp_idx);
bar(sorted_vals);
set(gca, 'XTickLabel', names(sorted_idx), 'TickLabelInterpreter', 'none');
ylabel('SampEn at Scale 1'); title('Scale 1 Entropy Ranking');
xtickangle(45); grid on;
bar(sorted_vals);
set(gca, 'XTickLabel', names(sorted_idx), 'TickLabelInterpreter', 'none');
ylabel('SampEn at Scale 1'); title('Scale 1 Entropy Ranking');
xtickangle(45); grid on;

% Scale behavior analysis
subplot(2,4,3);
scale_trends = zeros(num_sigs, 1);
for i = 1:num_sigs
    valid_scales = ~isnan(MSE_results(i,:));
    if sum(valid_scales) >= 5
        x_vals = scales(valid_scales);
        y_vals = MSE_results(i, valid_scales);
        p = polyfit(x_vals, y_vals, 1);
        scale_trends(i) = p(1); % slope
    else
        scale_trends(i) = NaN;
    end
end
bar(scale_trends);
set(gca, 'XTickLabel', names, 'TickLabelInterpreter', 'none');
ylabel('Scale Trend (slope)'); title('Scale Dependencies');
xtickangle(45); grid on;
yline(0, 'k--', 'LineWidth', 1);

% Computation times
subplot(2,4,4);
bar(computation_times);
set(gca, 'XTickLabel', names, 'TickLabelInterpreter', 'none');
ylabel('Time (s)'); title('Computation Times');
xtickangle(45); grid on;

% Parameter sensitivity - r values
subplot(2,4,5);
hold on;
r_fields = fieldnames(r_sensitivity);
for i = 1:length(r_fields)
    r_data = r_sensitivity.(r_fields{i});
    if isfield(r_data, 'white') && ~all(isnan(r_data.white))
        plot(1:length(r_data.white), r_data.white, 'o-', ...
            'DisplayName', sprintf('r=%.2f', config.r_values(i)));
    end
end
xlabel('Scale'); ylabel('SampEn'); title('r Parameter Sensitivity (White Noise)');
legend('Location', 'best'); grid on;

% Parameter sensitivity - m values  
subplot(2,4,6);
hold on;
m_fields = fieldnames(m_sensitivity);
for i = 1:length(m_fields)
    m_data = m_sensitivity.(m_fields{i});
    if isfield(m_data, 'white') && ~all(isnan(m_data.white))
        plot(1:length(m_data.white), m_data.white, 's-', ...
            'DisplayName', sprintf('m=%d', config.m_values(i)));
    end
end
xlabel('Scale'); ylabel('SampEn'); title('m Parameter Sensitivity (White Noise)');
legend('Location', 'best'); grid on;

% Edge case results
subplot(2,4,7);
edge_success_rates = [];
edge_labels = {};
for i = 1:numel(edge_names)
    edge_name = edge_names{i};
    edge_data = edge_results.(edge_name);
    
    if isstruct(edge_data)
        sub_names = fieldnames(edge_data);
        success = 0;
        for j = 1:numel(sub_names)
            if ~isnan(edge_data.(sub_names{j})(1))
                success = success + 1;
            end
        end
        success_rate = success / numel(sub_names);
    else
        success_rate = ~isnan(edge_data(1));
    end
    
    edge_success_rates(end+1) = success_rate * 100;
    edge_labels{end+1} = strrep(edge_name, '_', ' ');
end
bar(edge_success_rates);
set(gca, 'XTickLabel', edge_labels, 'TickLabelInterpreter', 'none');
ylabel('Success Rate (%)'); title('Edge Case Performance');
xtickangle(45); grid on;
ylim([0 105]);

% Signal complexity vs entropy at scale 1
subplot(2,4,8);
complexity_measure = zeros(num_sigs, 1);
for i = 1:num_sigs
    % Use lag-1 autocorrelation as complexity proxy
    complexity_measure(i) = 1 - abs(spectral_stats.(names{i}).lag1);
end
scatter(complexity_measure, scale1_values, 100, 'filled');
xlabel('Signal Complexity (1-|lag1|)'); ylabel('SampEn at Scale 1');
title('Complexity vs Entropy'); grid on;
for i = 1:num_sigs
    text(complexity_measure(i), scale1_values(i), names{i}, ...
         'FontSize', 8, 'HorizontalAlignment', 'center');
end

%% Theoretical Validation Analysis
fprintf('\n=== THEORETICAL VALIDATION RESULTS ===\n');
fprintf('%s\n', repmat('=', 1, 60));

% 1) Scale 1 ranking validation
fprintf('\n1. SCALE 1 ENTROPY RANKING\n');
fprintf('%s\n', repmat('-', 1, 30));
valid_idx = ~isnan(scale1_values);
valid_scale1 = scale1_values(valid_idx);
valid_names = names(valid_idx);
[~, observed_order] = sort(valid_scale1);
observed_names = valid_names(observed_order);

fprintf('Expected (low→high): %s\n', strjoin(expected_ranking, ' < '));
fprintf('Observed (low→high): %s\n', strjoin(observed_names, ' < '));

% Calculate ranking agreement
% Only compare signals that are present in both expected and observed rankings
valid_expected = [];
valid_observed = [];

for i = 1:length(expected_ranking)
    expected_signal = expected_ranking{i};
    % Find this signal in the observed ranking
    obs_idx = find(strcmpi(observed_names, expected_signal));
    if ~isempty(obs_idx)
        valid_expected(end+1) = i;  % Position in expected ranking
        valid_observed(end+1) = obs_idx; % Position in observed ranking
    end
end

if length(valid_expected) >= 3  % Need at least 3 signals for meaningful correlation
    ranking_correlation = corr(valid_expected', valid_observed', 'Type', 'Spearman');
else
    ranking_correlation = NaN; % Not enough signals for correlation
end
fprintf('Ranking correlation (Spearman ρ): %.3f\n', ranking_correlation);

% 2) Scale behavior validation
fprintf('\n2. SCALE BEHAVIOR VALIDATION\n');
fprintf('%s\n', repmat('-', 1, 30));
behavior_checks_passed = 0;
total_behavior_checks = 0;

for i = 1:num_sigs
    signal_name = names{i};
    if isKey(scale_behaviors, signal_name)
        expected_behavior = scale_behaviors(signal_name);
        total_behavior_checks = total_behavior_checks + 1;
        
        % Analyze actual behavior
        valid_scales = ~isnan(MSE_results(i,:));
        if sum(valid_scales) >= 5
            trend_slope = scale_trends(i);
            
            behavior_match = false;
            switch expected_behavior
                case 'decreasing'
                    behavior_match = trend_slope < -0.01;
                case 'increasing'  
                    behavior_match = trend_slope > 0.01;
                case 'flat_low'
                    behavior_match = abs(trend_slope) < 0.01 && mean(MSE_results(i,:), 'omitnan') < 0.5;
                case 'mild_decrease'
                    behavior_match = trend_slope < 0 && trend_slope > -0.05;
                case 'hump'
                    % Check for peak somewhere in middle scales
                    [~, peak_idx] = max(MSE_results(i,:));
                    behavior_match = peak_idx > 3 && peak_idx < num_scales-3;
            end
            
            if behavior_match
                behavior_checks_passed = behavior_checks_passed + 1;
                status = '✓';
            else
                status = '✗';
            end
            
            fprintf('%s %-10s: expected %-12s, slope=%.4f %s\n', ...
                status, signal_name, expected_behavior, trend_slope, status);
        else
            fprintf('? %-10s: insufficient valid data\n', signal_name);
        end
    end
end

behavior_pass_rate = behavior_checks_passed / total_behavior_checks * 100;
fprintf('\nScale behavior validation: %d/%d passed (%.0f%%)\n', ...
        behavior_checks_passed, total_behavior_checks, behavior_pass_rate);

% 3) Parameter sensitivity validation
fprintf('\n3. PARAMETER SENSITIVITY\n');
fprintf('%s\n', repmat('-', 1, 30));

% Check r parameter sensitivity
r_sensitivity_ok = true;
fprintf('r parameter sensitivity:\n');
for i = 1:length(config.r_values)-1
    r1_key = sprintf('r_%03d', round(config.r_values(i)*100));
    r2_key = sprintf('r_%03d', round(config.r_values(i+1)*100));
    
    if isfield(r_sensitivity, r1_key) && isfield(r_sensitivity, r2_key)
        if isfield(r_sensitivity.(r1_key), 'white') && isfield(r_sensitivity.(r2_key), 'white')
            vals1 = r_sensitivity.(r1_key).white;
            vals2 = r_sensitivity.(r2_key).white;
            
            if ~all(isnan(vals1)) && ~all(isnan(vals2))
                mean_diff = mean(vals2 - vals1, 'omitnan');
                fprintf('  r=%.2f→%.2f: Δ=%.3f %s\n', ...
                    config.r_values(i), config.r_values(i+1), mean_diff, ...
                    tern(abs(mean_diff) > 0.01, '(sensitive)', '(stable)'));
            end
        end
    end
end

% Check m parameter sensitivity  
fprintf('m parameter sensitivity:\n');
m_fields = fieldnames(m_sensitivity);
if length(m_fields) >= 2
    for i = 1:length(m_fields)-1
        m_data1 = m_sensitivity.(m_fields{i});
        m_data2 = m_sensitivity.(m_fields{i+1});
        
        if isfield(m_data1, 'white') && isfield(m_data2, 'white')
            vals1 = m_data1.white;
            vals2 = m_data2.white;
            
            if ~all(isnan(vals1)) && ~all(isnan(vals2))
                mean_diff = mean(vals2 - vals1, 'omitnan');
                fprintf('  m=%d→%d: Δ=%.3f %s\n', ...
                    config.m_values(i), config.m_values(i+1), mean_diff, ...
                    tern(abs(mean_diff) > 0.05, '(sensitive)', '(stable)'));
            end
        end
    end
end

% 4) Edge case robustness
fprintf('\n4. EDGE CASE ROBUSTNESS\n');
fprintf('%s\n', repmat('-', 1, 30));
edge_pass_count = 0;
for i = 1:length(edge_success_rates)
    if edge_success_rates(i) >= 100
        edge_pass_count = edge_pass_count + 1;
        status = '✓';
    else
        status = '✗';
    end
    fprintf('%s %-20s: %.0f%% success\n', status, edge_labels{i}, edge_success_rates(i));
end
edge_robustness = edge_pass_count / length(edge_success_rates) * 100;

%% Final Assessment
fprintf('\n=== OVERALL VALIDATION ASSESSMENT ===\n');
fprintf('%s\n', repmat('=', 1, 60));

% Compute overall validation score
ranking_score = max(0, min(100, (ranking_correlation + 1) * 50)); % Convert -1:1 to 0:100
behavior_score = behavior_pass_rate;
robustness_score = edge_robustness;
overall_score = (ranking_score + behavior_score + robustness_score) / 3;

fprintf('\nValidation Component Scores:\n');
fprintf('• Signal ranking correlation: %.0f%% (ρ=%.3f)\n', ranking_score, ranking_correlation);
fprintf('• Scale behavior validation: %.0f%%\n', behavior_score);
fprintf('• Edge case robustness: %.0f%%\n', robustness_score);
fprintf('• Overall validation score: %.0f%%\n\n', overall_score);

% Performance assessment
mean_computation_time = mean(computation_times, 'omitnan');
total_computation_time = sum(computation_times);
fprintf('Performance Metrics:\n');
fprintf('• Average computation time: %.2fs per signal\n', mean_computation_time);
fprintf('• Total analysis time: %.1fs\n', total_computation_time);

% Final recommendation
if overall_score >= 90
    fprintf('\n✓ EXCELLENT: Algorithm passes comprehensive validation\n');
    fprintf('  Recommendation: Ready for production use\n');
elseif overall_score >= 75
    fprintf('\n✓ GOOD: Algorithm shows strong performance\n');
    fprintf('  Recommendation: Suitable for most applications\n');
elseif overall_score >= 60
    fprintf('\n⚠ MODERATE: Algorithm has some limitations\n');
    fprintf('  Recommendation: Use with awareness of specific limitations\n');
else
    fprintf('\n✗ POOR: Algorithm fails multiple validation criteria\n');
    fprintf('  Recommendation: Requires improvement before use\n');
end

% Specific insights and recommendations
fprintf('\nKey Insights:\n');
if ranking_correlation > 0.7
    fprintf('• Strong correlation with theoretical entropy ranking\n');
else
    fprintf('• Weak correlation with theoretical expectations - investigate discrepancies\n');
end

if behavior_pass_rate > 80
    fprintf('• Scale behaviors match theoretical expectations well\n');
else
    fprintf('• Several signals show unexpected scale behaviors - review implementation\n');
end

if mean_computation_time < 1.0
    fprintf('• Excellent computational efficiency\n');
elseif mean_computation_time < 5.0
    fprintf('• Good computational efficiency\n');
else
    fprintf('• Slow computation - consider optimization\n');
end

% Signal-specific warnings
problematic_signals = [];
for i = 1:num_sigs
    if all(isnan(MSE_results(i,:)))
        problematic_signals{end+1} = names{i};
    end
end

if ~isempty(problematic_signals)
    fprintf('• WARNING: Failed to compute MSE for: %s\n', strjoin(problematic_signals, ', '));
end

fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('Validation completed successfully!\n');
fprintf('Review plots and results for detailed analysis.\n');

%% Helper Functions

function y = tern(condition, true_val, false_val)
    if condition
        y = true_val;
    else
        y = false_val;
    end
end