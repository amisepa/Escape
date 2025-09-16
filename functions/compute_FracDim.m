function [D, SD, info] = compute_FracDim(data, varargin)
% Fractal dimension via box counting with outlier-aware options.
%
%   [D, SD, info] = compute_FracDim(data, ...
%       'RejectBursts', true, 'WinFrac', 0.02, 'ZThresh', 6, ...
%       'RobustFit', 'theilsen', 'ScaleTrimIQR', true, ...
%       'MinBoxesPerCol', 1, 'MinJ', [], 'MaxJ', [], ...
%       'Parallel', true, 'Progress', true, 'MinScales', 6)
%
% Inputs
%   data            : [n_channels x n_samples] (vector → 1 channel)
%   Name-Value pairs:
%     'RejectBursts'   : true|false, MAD z-score burst rejection (default true)
%     'WinFrac'        : window length for burst detection as fraction of x-axis (default 0.02)
%     'ZThresh'        : MAD z threshold (default 6)
%     'RobustFit'      : 'theilsen' | 'ols'  (default 'theilsen')
%                        Note: 'huber' is mapped to 'theilsen' for portability.
%     'ScaleTrimIQR'   : true|false, trim scale-wise slope outliers by IQR (default true)
%     'MinBoxesPerCol' : minimum vertical boxes per x-column (default 1)
%     'MinJ'           : minimum dyadic level (default [], auto)
%     'MaxJ'           : maximum dyadic level (default [], auto)
%     'Parallel'       : true|false, parfor over channels (default true)
%     'Progress'       : true|false, show progress (default true)
%                        • Parallel==true  & Progress==true → text-only (parfor-safe)
%                        • Parallel==false & Progress==true → text + waitbar (fallback to text if headless)
%     'MinScales'      : minimum #valid scales for robust fit (default 6)
%
% Outputs
%   D   : fractal dimension per channel [n_channels x 1]
%   SD  : standard error of slope per channel [n_channels x 1]
%   info: diagnostics (fields per-channel: .r, .n, .validScales, .usedX, .usedY,
%         .cleanMask, .Jvec, .fitMethod, .nUsedScales, .flags; plus .params)
%
% -------------------------------------------------------------------------
% Copyright (C) 2025
% EEGLAB Escape Plugin — Author: Cedric Cannard
% License: GNU GPL v2 or later
% -------------------------------------------------------------------------

% ---------- Parse inputs ----------
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
p.addParameter('Parallel', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Progress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MinScales', 6, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.parse(data, varargin{:});
opts = p.Results;

% Map 'huber' to 'theilsen' for portability
if strcmpi(opts.RobustFit,'huber')
    warning('compute_FracDim:HuberNotImplemented', ...
        'Huber robust fit not provided in this portable build. Using Theil–Sen instead.');
    opts.RobustFit = 'theilsen';
end

% ---------- Shape ----------
if size(data,1) > size(data,2)
    data = data.';  % -> [n_channels x n_samples]
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
info.fitMethod   = cell(nchan,1);
info.rejectBursts= repmat(logical(opts.RejectBursts), nchan, 1);
info.nUsedScales = zeros(nchan,1);
info.flags       = cell(nchan,1);
info.params      = opts;

% ---------- Progress header ----------
if opts.Progress
    if opts.Parallel
        fprintf('FracDim: %d channel(s) | robust=%s | parallel=on (text only)\n', ...
            nchan, lower(string(opts.RobustFit)));
        fprintf('Progress:\n');
    else
        fprintf('FracDim: %d channel(s) | robust=%s | parallel=off (text + waitbar)\n', ...
            nchan, lower(string(opts.RobustFit)));
    end
end

% ---------- Iterate channels ----------
if opts.Parallel && ~isempty(ver('parallel'))
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

        if opts.Progress
            fprintf('  ch %3d/%3d: D=%7.5f  SE=%7.5f\n', ch, nchan, D(ch), SD(ch));
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
    useWB = opts.Progress && usejava('desktop');
    hWB = [];
    if useWB
        try hWB = waitbar(0,'Computing Fractal Dimension...','Name','compute_FracDim'); catch, hWB = []; end
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

        if opts.Progress
            fprintf('  ch %3d/%3d: D=%7.5f  SE=%7.5f\n', ch, nchan, D(ch), SD(ch));
            if ~isempty(hWB) && isvalid(hWB)
                try, waitbar(ch/nchan, hWB, sprintf('Computing Fractal Dimension... (%d/%d)', ch, nchan)); end
            end
        end
    end
    if ~isempty(hWB) && isvalid(hWB), try, close(hWB); end, end
end
end

%% ===================== Single-channel worker ===================== %%
function [dimension, standard_dev, out] = compute_fractal_volatility_single(sig, opts)
sig = sig(:).';           % row
N   = numel(sig);

% Build [x,y] in unit square
x = (1:N).';
y = sig(:);

% Remove NaNs/Infs early
bad = ~isfinite(y);
if any(bad), y(bad) = []; x(bad) = []; N = numel(y); end

% Normalize to unit square
x = (x - min(x)) / max(eps, (max(x)-min(x)));
y = (y - min(y)) / max(eps, (max(y)-min(y)));

% Guard constant signals
if isempty(y) || max(y) == min(y)
    dimension = NaN; standard_dev = NaN;
    out = struct('r',[],'n',[],'validScales',[],'usedX',[],'usedY',[], ...
                 'cleanMask',~bad, 'Jvec',[], 'fitMethod','none', ...
                 'nUsedScales',0,'flags',struct('fewScales',true,'highSE',false,'constant',true));
    return
end

% Optional burst rejection via sliding MAD z-score
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
    xg = x; yg = y; maskGood = true(N,1);  % fallback if too much removed
end

% Dyadic scale selection
dx = diff(xg); dx(dx<=0) = [];
if isempty(dx), dx = 1/N; end
minwidth = floor(abs(log2(max(eps, min(dx))))) - 1;
minwidth = max(minwidth, 2);
Jmax_auto = minwidth; Jmin_auto = 1;
Jmin = Jmin_auto; if ~isempty(opts.MinJ), Jmin = max(1, opts.MinJ); end
Jmax = Jmax_auto; if ~isempty(opts.MaxJ), Jmax = min(opts.MaxJ, Jmax_auto); end
Jvec = Jmin:Jmax;

% Box counting per scale
n = zeros(numel(Jvec),1);
r = zeros(numel(Jvec),1);
for ii = 1:numel(Jvec)
    j = Jvec(ii);
    width = 2^(-j);
    n(ii) = count_boxes_in_columns(xg, yg, width, opts.MinBoxesPerCol);
    r(ii) = width;
end

% Valid scales
valid = n > 0 & r > 0 & isfinite(n) & isfinite(r);
rv = r(valid); nv = n(valid);
if numel(rv) < 2
    dimension = NaN; standard_dev = NaN;
    out = struct('r',r,'n',n,'validScales',valid,'usedX',[],'usedY',[], ...
                 'cleanMask',maskGood,'Jvec',Jvec,'fitMethod','none', ...
                 'nUsedScales',0,'flags',struct('fewScales',true,'highSE',false,'constant',false));
    return
end

% Optional scale-trim via IQR on local slopes
if opts.ScaleTrimIQR && numel(rv) >= 5
    s = -gradient(log(nv)) ./ gradient(log(rv));
    IQRs = iqr(s);
    keep = abs(s - median(s)) <= IQRs/2;
    if numel(keep) == numel(rv)
        rv = rv(keep); nv = nv(keep);
    end
end

x2 = log(rv); y2 = log(nv);

% Standardize for robust fit stability
mx = mean(x2); sx = std(x2); if sx==0, sx = 1; end
my = mean(y2); sy = std(y2); if sy==0, sy = 1; end
xz = (x2 - mx) / sx;  xz = xz(:);
yz = (y2 - my) / sy;  yz = yz(:);

% Choose fit method
forceOLS = strcmpi(opts.RobustFit,'ols');
smallN   = numel(xz) < opts.MinScales || rank([ones(numel(xz),1) xz]) < 2;

if forceOLS || smallN
    % OLS on (x2,y2)
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
    % Theil–Sen on standardized vars; SE from OLS on standardized residuals
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

%% ===================== Helpers ===================== %%
function boxcount = count_boxes_in_columns(xg, yg, width, minBoxes)
% Count vertical boxes needed inside each x column of given width.
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
K = N*(N-1)/2; %#ok<NASGU> 
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
