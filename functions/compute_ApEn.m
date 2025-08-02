% Computes the approximate entropy (AE)
%
% Inputs:
%   x   - univariate signal - a vector of size 1 x n (the number of sample points)
%   m   - embedding dimension (default = 2)
%   r   - threshold/tolerance (.15 of the signal' sd)
% 
% Outputs:
%   ae - approximate entropy
%
% Cedric Cannard


function ae = compute_ApEn(signal, m, r)

% Defaults
if ~exist('m', 'var'), m = 2; end
if ~exist('r', 'var'), r = .15*std(signal); end
if ~exist('tau', 'var'), tau = 1; end

N = length(signal);
result = zeros(1,2);

for j = 1:2
    m = m+j-1;
    phi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);
    
    % setting up data matrix
    parfor i = 1:m
        dataMat(i,:) = signal(i:N-m+i);
    end
    
    % counting similar patterns using distance calculation
    parfor i = 1:N-m+1
        
        if sum( isnan( dataMat(:,i) ) ) == 0
        
            tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));
            boolMat = any( (tempMat > r),1);
            phi(i) = sum(~boolMat)/(N-m+1);
        else

	        %discarding blocks with nan values
            phi(i) = nan;
        end
    end
    
    % summing over the counts
    result(j) = log(mean(phi, 'omitnan'));    
end

ae = result(1) - result(2) ;

end