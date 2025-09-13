function [D, SD, info] = compute_FracDim(data, varargin)
% Fractal dimension via box counting with outlier-aware options.
%
%   [D, SD, info] = compute_FracDim(data, 'Parallel', true, 'RejectBursts', true, ...)
%
%   Inputs
%     - data          : matrix [n_channels x n_samples]. If a vector, treated as 1 channel.
%     Name-Value pairs:
%       'RejectBursts'   : true|false, sliding MAD z-score burst rejection (default true)
%       'WinFrac'        : window length for burst detection as fraction of [0,1] x-axis (default 0.02)
%       'ZThresh'        : MAD z threshold to flag bursts (default 6)
%       'RobustFit'      : 'theilsen' | 'ols' (default 'theilsen')
%                          Note: if 'huber' is provided, it is mapped to 'theilsen' for portability.
%       'ScaleTrimIQR'   : true|false, trim scale-wise slope outliers by IQR rule (default true)
%       'MinBoxesPerCol' : minimum vertical boxes to accept per column (default 1)
%       'MinJ'           : minimum dyadic level to consider, [] for auto (default [])
%       'MaxJ'           : maximum dyadic level to consider, [] for auto (default [])
%       'Parallel'       : true|false, enable parfor over channels (default false)
%       'Verbose'        : true|false, print per-channel info and progress (default true)
%       'MinScales'      : minimum number of valid scales required to use the robust fit (default 6)
%
%   Outputs
%     - D   : fractal dimension per channel [n_channels x 1]
%     - SD  : standard error of slope per channel [n_channels x 1]
%     - info: struct with diagnostics
%         .r, .n, .validScales, .usedX, .usedY, .cleanMask, .Jvec per channel (cell arrays)
%         .fitMethod, .rejectBursts, .nUsedScales, .flags, .params
%
%   Rationale and references
%     Short, high-power bursts can distort scale-based estimates even if they occupy a
%     tiny fraction of the recording. Pre-removal of such outliers and robust fitting of
%     scaling relationships stabilizes the slope estimate. Methodological guidance taken
%     from wavelet-leader multifractal outlier handling and robust estimation ideas:
%       Dumeur, M., Palva, J. M., & Ciuciu, P. (2025). Outlier detection and removal in
%       multifractal analysis of electrophysiological brain signals. EURASIP Journal on
%       Advances in Signal Processing, 2025(1), 35.
%
%   Copyright (C) Cedric Cannard 2025
%   Escape EEGLAB Plugin  https://github.com/amisepa/Escape
%
%   This code adapts general recommendations about detecting and downweighting outlier
%   segments before scale-exponent fitting. It does not implement the paper's full
%   p-leader segmentation, but mirrors its intent for box counting.

% Parse inputs
p = inputParser;
p.addRequired('data', @(x) isnumeric(x) && ismatrix(x));
p.addParameter('RejectBursts', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('WinFrac',      0.02, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ZThresh',      6,    @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('RobustFit',    'theilsen', @(s) ischar(s) || isstring(s));
p.addParameter('ScaleTrimIQR', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MinBoxesPerCol', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('MinJ', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('MaxJ', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
p.addParameter('Parallel', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Verbose',  true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('MinScales', 6, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.parse(data, varargin{:});
opts = p.Results;

% Map 'huber' to 'theilsen' for portability
if strcmpi(opts.RobustFit,'huber')
    warning('compute_FracDim:HuberNotImplemented', ...
        'Huber robust fit not provided in this portable build. Using Theil-Sen instead.');
    opts.RobustFit = 'theilsen';
end

% Ensure [n_channels x n_samples]
if size(data,1) > size(data,2)
    data = data.';  % transpose to [n_channels x n_samples]
end
[nchan, ~] = size(data);
D  = nan(nchan,1);
SD = nan(nchan,1);

% Pre-alloc info containers
info = struct();
info.r           = cell(nchan,1);
info.n           = cell(nchan,1);
info.validScales = cell(nchan,1);
info.usedX       = cell(nchan,1);
info.usedY       = cell(nchan,1);
info.cleanMask   = cell(nchan,1);
info.Jvec        = cell(nchan,1);
info.fitMethod   = cell(nchan,1);  % per-channel actual method used
info.rejectBursts= repmat(logical(opts.RejectBursts), nchan, 1);
info.nUsedScales = zeros(nchan,1);
info.flags       = cell(nchan,1);
info.params      = opts;

if opts.Verbose
    fprintf('Fractal volatility via box counting on %d channels...\n', nchan);
end

% Iterate channels, parallel or serial
if opts.Parallel
    % Collectors for parfor
    r_all     = cell(nchan,1);
    n_all     = cell(nchan,1);
    valid_all = cell(nchan,1);
    x_all     = cell(nchan,1);
    y_all     = cell(nchan,1);
    mask_all  = cell(nchan,1);
    J_all     = cell(nchan,1);
    method_all= cell(nchan,1);
    nscale_all= zeros(nchan,1);
    flags_all = cell(nchan,1);

    parfor ch = 1:nchan %#ok<PFOUS>
        [D(ch), SD(ch), chInfo] = compute_fractal_volatility_single(data(ch,:), opts);
        r_all{ch}      = chInfo.r;
        n_all{ch}      = chInfo.n;
        valid_all{ch}  = chInfo.validScales;
        x_all{ch}      = chInfo.usedX;
        y_all{ch}      = chInfo.usedY;
        mask_all{ch}   = chInfo.cleanMask;
        J_all{ch}      = chInfo.Jvec;
        method_all{ch} = chInfo.fitMethod;
        nscale_all(ch) = chInfo.nUsedScales;
        flags_all{ch}  = chInfo.flags;

        if opts.Verbose
            fprintf('Channel %3d/%3d  D=%6.4f  SE=%6.4f\n', ch, nchan, D(ch), SD(ch));
        end
    end

    % Repackage
    info.r            = r_all;
    info.n            = n_all;
    info.validScales  = valid_all;
    info.usedX        = x_all;
    info.usedY        = y_all;
    info.cleanMask    = mask_all;
    info.Jvec         = J_all;
    info.fitMethod    = method_all;
    info.nUsedScales  = nscale_all;
    info.flags        = flags_all;

else
    if exist('progressbar','file')
        progressbar('computing fractal dimension/volatility (serial)')
    end
    for ch = 1:nchan
        [D(ch), SD(ch), chInfo] = compute_fractal_volatility_single(data(ch,:), opts);
        info.r{ch}           = chInfo.r;
        info.n{ch}           = chInfo.n;
        info.validScales{ch} = chInfo.validScales;
        info.usedX{ch}       = chInfo.usedX;
        info.usedY{ch}       = chInfo.usedY;
        info.cleanMask{ch}   = chInfo.cleanMask;
        info.Jvec{ch}        = chInfo.Jvec;
        info.fitMethod{ch}   = chInfo.fitMethod;
        info.nUsedScales(ch) = chInfo.nUsedScales;
        info.flags{ch}       = chInfo.flags;

        if opts.Verbose
            fprintf('Channel %3d/%3d  D=%6.4f  SE=%6.4f\n', ch, nchan, D(ch), SD(ch));
        end
        if exist('progressbar','file')
            progressbar(ch / nchan)
        end
    end
end
end


%% Single-channel worker
function [dimension, standard_dev, out] = compute_fractal_volatility_single(sig, opts)
sig = sig(:).';           % row
N   = numel(sig);

% Build [x,y] embedded in unit square
x = (1:N).';
y = sig(:);

% Remove NaNs/Infs early
bad = ~isfinite(y);
if any(bad), y(bad) = []; x(bad) = []; N = numel(y); end

% Normalize to unit square
x = (x - min(x)) / max(eps, (max(x)-min(x)));
y = (y - min(y)) / max(eps, (max(y)-min(y)));

% Guard constant signals (no vertical extent)
if isempty(y) || max(y) == min(y)
    dimension = NaN; standard_dev = NaN;
    out = struct('r',[],'n',[],'validScales',[],'usedX',[],'usedY',[], ...
                 'cleanMask',~bad, 'Jvec',[], 'fitMethod','none', ...
                 'nUsedScales',0,'flags',struct('fewScales',true,'highSE',false,'constant',true));
    return
end

% Optional burst rejection using sliding MAD z-score on y
maskGood = true(N,1);
if opts.RejectBursts
    win = max(5, round(opts.WinFrac * N));
    y_med = movmedian(y, win, 'Endpoints','shrink');
    y_mad = movmad(y, win, 1, 'Endpoints','shrink');
    mad_eps = median(y_mad(y_mad>0));
    if isempty(mad_eps) || mad_eps==0, mad_eps = 1e-6; end
    z = abs(y - y_med) ./ max(y_mad, mad_eps);
    maskGood = z <= opts.ZThresh;
end

xg = x(maskGood); yg = y(maskGood);
[xg, idx] = sort(xg, 'ascend'); yg = yg(idx);
if numel(xg) < 10
    % fall back to raw if cleaning removed too much
    xg = x; yg = y; maskGood = true(N,1);
end

% Determine dyadic scales based on min spacing of x
dx = diff(xg); dx(dx<=0) = [];
if isempty(dx), dx = 1/N; end
minwidth = floor(abs(log2(max(eps, min(dx))))) - 1;
minwidth = max(minwidth, 2);
Jmax_auto = minwidth; Jmin_auto = 1;
Jmin = Jmin_auto; if ~isempty(opts.MinJ), Jmin = max(1, opts.MinJ); end
Jmax = Jmax_auto; if ~isempty(opts.MaxJ), Jmax = min(opts.MaxJ, Jmax_auto); end
Jvec = Jmin:Jmax;

% Box counting over columns for each scale
n = zeros(numel(Jvec),1);
r = zeros(numel(Jvec),1);
for ii = 1:numel(Jvec)
    j = Jvec(ii);
    width = 2^(-j);
    n(ii) = count_boxes_in_columns(xg, yg, width, opts.MinBoxesPerCol);
    r(ii) = width;
end

% Assemble log arrays and drop degenerate scales
valid = n > 0 & r > 0 & isfinite(n) & isfinite(r);
rv = r(valid); nv = n(valid);
if numel(rv) < 2
    dimension = NaN; standard_dev = NaN;
    out = struct('r',r,'n',n,'validScales',valid,'usedX',[],'usedY',[], ...
                 'cleanMask',maskGood,'Jvec',Jvec,'fitMethod','none', ...
                 'nUsedScales',0,'flags',struct('fewScales',true,'highSE',false,'constant',false));
    return
end

% Optional trimming of scale-wise slopes by IQR
if opts.ScaleTrimIQR && numel(rv) >= 5
    s = -gradient(log(nv)) ./ gradient(log(rv));
    IQRs = iqr(s);
    keep = abs(s - median(s)) <= IQRs/2;
    if numel(keep) == numel(rv)
        rv = rv(keep); nv = nv(keep);
    end
end

x2 = log(rv); y2 = log(nv);

% Standardize and ensure column vectors
mx = mean(x2); sx = std(x2); if sx==0, sx = 1; end
my = mean(y2); sy = std(y2); if sy==0, sy = 1; end
xz = (x2 - mx) / sx;  xz = xz(:);
yz = (y2 - my) / sy;  yz = yz(:);

% Decide method
forceOLS = strcmpi(opts.RobustFit,'ols');
smallN   = numel(xz) < opts.MinScales || rank([ones(numel(xz),1) xz]) < 2;

if forceOLS || smallN
    % OLS on original variables
    X = [ones(size(x2)) x2];
    beta = X \ y2;
    yhat = X*beta;
    e = y2 - yhat;
    C = pinv(X.'*X);
    sigma2 = (e.'*e) / max(1, numel(y2)-2);
    se = sqrt(sigma2) .* sqrt(diag(C));
    slope    = beta(2);
    slope_se = se(2);
    methodUsed = 'ols';
else
    % Theil–Sen on standardized variables, then back-transform
    m = theil_sen(xz, yz);
    b = median(yz - m*xz);
    slope    = m * sy / sx;
    % SE via OLS on standardized data
    Xz = [ones(size(xz)) xz];
    ez = yz - (b + m*xz);
    sigma2z = (ez'*ez) / max(1, numel(yz)-2);
    Cz = pinv(Xz.'*Xz);
    se_z = sqrt(sigma2z) * sqrt(diag(Cz));
    slope_se = se_z(2) * sy / sx;
    methodUsed = 'theilsen';
end

dimension    = -slope;
standard_dev = slope_se;

% Pack diagnostics
out = struct();
out.r = r;
out.n = n;
out.validScales = valid;
out.usedX = x2;
out.usedY = y2;
out.cleanMask = maskGood;
out.Jvec = Jvec;
out.fitMethod = methodUsed;
out.nUsedScales = numel(x2);
out.flags = struct( ...
    'fewScales', numel(x2) < opts.MinScales, ...
    'highSE',    ~isnan(standard_dev) && standard_dev > 0.05, ...
    'constant',  false);
end

%% Helpers
function boxcount = count_boxes_in_columns(xg, yg, width, minBoxes)
% Count vertical boxes needed inside each x column of width.
nCols = max(1, 2^round(-log2(width)));
col = min(floor(xg / width) + 1, nCols);
boxcount = 0;
for c = 1:nCols
    idx = (col == c);
    if ~any(idx), continue; end
    yc = yg(idx);
    if numel(yc) == 1
        boxcount = boxcount + 1;
    else
        rawcount = (max(yc) - min(yc)) / width;
        rawcount = rawcount + rem(min(yc), width);
        cnt = ceil(rawcount);
        if cnt >= minBoxes
            boxcount = boxcount + cnt;
        end
    end
end
end

function m = theil_sen(x, y)
x = x(:); y = y(:);
N = numel(x);
if N < 2, m = 0; return, end
K = N*(N-1)/2;
sl = zeros(K,1);
t = 0;
for i = 1:N-1
    dx = x(i+1:end) - x(i);
    dy = y(i+1:end) - y(i);
    k = dx ~= 0 & isfinite(dx) & isfinite(dy);
    nn = sum(k);
    if nn > 0
        sl(t+1:t+nn) = dy(k) ./ dx(k);
        t = t + nn;
    end
end
if t == 0
    m = (y(end)-y(1)) / max(eps, (x(end)-x(1)));
else
    m = median(sl(1:t));
end
end
