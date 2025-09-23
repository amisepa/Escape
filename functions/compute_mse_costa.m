function [ e, A, B ] = compute_MSE_costa( x, m, r, tau, coarsing)
% MULTISCALE SAMPLE ENTROPY (MSE)
%
% Based on "Multiscale entropy analysis of biological signals"
% By Madalena Costa, Ary L. Goldberger, and C.-K. Peng
% Published on 18 February 2005 in Phys. Rev. E 71, 021906.
%
% Author: Cedric Cannard

switch nargin
    case 1
        m = 2;
        r = 0.15;
        tau = 1;
        coarsing = 'std';
    case 2
        r = 0.15;
        tau = 1;
    case 3
        tau = 1;
end

% skip 1st scale for std and var
if tau == 1
    e = NaN;
    A = NaN;
    B = NaN;
    return
end

% coarse signal
if strcmpi(coarsing,'mean')
    y = mean(buffer(x(:), tau), 1, 'omitmissing');
elseif strcmpi(coarsing,'median')
    y = median(buffer(x(:), tau), 1, 'omitmissing');
elseif strcmpi(coarsing,'std')
    y = std(buffer(x(:), tau), 1, 'omitmissing');
elseif strcmpi(coarsing,'var')
    y = variance(buffer(x(:), tau), 1, 'omitmissing');
else
    error("coarsing method unrecognized. Must be 'mean', 'median', 'std', or 'var'.")
end

% (m+1)-element sequences
X = buffer(y, m + 1, m, 'nodelay')';

% matching (m+1)-element sequences
A = sum(pdist(X, 'chebychev') < r * std(x, 1,'omitmissing'));

% matching m-element sequences
X = X(:, 1:m);
B = sum(pdist(X, 'chebychev') < r * std(x, 1,'omitmissing'));

% take log
if A == 0 || B == 0
    e = NaN;
    return
end
e = log(B / A);

end

