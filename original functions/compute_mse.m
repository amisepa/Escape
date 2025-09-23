%% Computes multiscale entropy (MSE)
%
% INPUTS:
%   signal: univariate signal - a vector of size 1 x N (the number of sample points)
%   m: embedding dimension
%   r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - 
%       because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
%   tau: time lag (it is usually equal to 1)
%   coarseType: 'Mean', 'Standard deviation' (default), 'Variance'
%   nScales: the number of scale factors (default = 15)
%   filtData: bandpass filter each scale factor to control for spectral 
%       bias (1) or not (0; default)
%   fs: sample rate (Hz)
%
% OUTPUTS:
%   mse: entropy values for each scale factor
%   scales: lower and upper frequency bounds of each time scale
%
% EXAMPLE:
%   [mse, scales] = compute_mse(signal, m, r, tau, coarseType, nScales, filtData, fs)
% 
% Ref:
%   [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on
%   Standard Deviation for Biomedical Signal Analysis", Medical & Biological
%   Engineering & Computing, 2016.
%   [2] M. Costa, A. Goldberger, and C.K. Peng, "Multiscale Entropy Analysis
%   of Complex Physiologic Time Series", Physical review letters, vol. 89, no. 6, p. 068102, 2002.
%
% Cedric Cannard, 2022

function [mse, scales] = compute_mse(signal, m, r, tau, coarseType, nScales, filtData, fs)

% Max scale factor cannot be greater than Nyquist frequency
nf = fs/2;
if nScales >= nf
    warning(["Scale factor cannot be as high as the Nyquist frequency. Lowering it to " num2str(nf-1) ]);
end

% Signal is centered and normalized to have SD = 1
signal = signal-mean(signal);
signal = signal./std(signal);

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
end

mse = nan(1,nScales);
scales = nan(2,nScales);
parfor iScale = 1:nScales
    
    % Make copy of signal in case it is bandpass-filtered at each scale
    sig = signal;
    
    % Scale factor bounds
    upperBound = (1/iScale).*nf + .05*((1./iScale).*nf);
    lowerBound = (1/(iScale+1)).*nf - .05*((1./(iScale+1)).*nf);

    % Display the corresponding frequencies for this scale
    scales(:,iScale) = [round(lowerBound,3) round(upperBound,3) ];
%     disp(['scale ' num2str(iScale) ': ' num2str(round(lowerBound,3)) ' - ' num2str(round(upperBound,3)) ' Hz']);
    fprintf('   scale %d \n', iScale)

    % Bandpass filter outside these bounds to control for spectral bias (if selected)
    if filtData
        disp('Applying bandpasss-filter remove spectral bias from outside these bounds.')

        if iScale > 1
            if upperBound-lowerBound > .05*nf    % for broad passbands: Chebyshev Type I zero-phase bandpass
                [b,a] = cheby1(4,1,upperBound/nf,'low');
                % figure; freqz(b,a);
                sig = filtfilt(b,a,sig);
                [b,a] = cheby1(4,1,lowerBound/nf,'high');
                sig = filtfilt(b,a,sig);
            else                                 % Butterworth zero-phase bandpass
                [b,a] = butter(10,upperBound/nf,'low');
                sig = filtfilt(b,a,sig);
                [b,a] = butter(10,lowerBound/nf,'high');
                sig = filtfilt(b,a,sig);
            end
        else
            [b,a]= butter(10,lowerBound/nf,'high');   % only Butterworth highpass for scale 1
            % upperBound = nf;
            sig = filtfilt(b,a,sig);
        end

%         % Visualize filter effect on the power spectrum
%         [psd,f] = pwelch(sig,[],[],[],fs);
%         plot(f,psd); hold on;
%         title(num2str(iScale)); legend([num2str(lowerBound) '-' num2str(upperBound)]);

    end

    y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);

    switch coarseType
        case 'Mean'
            sig = mean(y,'omitnan');
        case 'SD'
            sig = std(y,'omitnan');
        case 'Variance'
            sig = var(y,'omitnan');
    end

    mse(:,iScale) = compute_se(sig, m, r, tau);

end

% % Remove NaN scales
% idx = isnan(mse);
% mse(idx) = [];
% scales(idx) = [];
