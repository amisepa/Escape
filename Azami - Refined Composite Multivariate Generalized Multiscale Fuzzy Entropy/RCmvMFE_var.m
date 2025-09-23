function Out_mvMFE=RCmvMFE_var(X,sm,sr,n,stau,Scale)
%
% To improve the stabily and reliabilty of multivaraite multiscale fuzzy entropy (mvMFE) based on variance (mvMFE_var), especially for short signals, we proposed refined composite mvMFE_var (RCmvMFE_var). In RCmvMFE_var, for scale
%factor tau, tau different multivaraite time series, corresponding to different starting points of the coarse-graining process are created and the RCmvMFE_var value is defined based on
%the averages of the total number of m- and m+1- dimensional matched vector pairs of those shifted sequences.
%
%
% Inputs:
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% sm: scalar embedding value
% sr: scalar threshold value (it is usually equal to 0.15)
% n: fuzzy power (it is usually equal to 2)
% stau: scalar time lag  value (it is usually equal to 1)
% Scale: the number of scale factors
%
% Output:
% Out_mvMFE: a scalar quantity
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16



%%
% Because multi-channel signals may have different amplitude ranges, the distances calculated on embedded vectors may be biased
%toward the largest amplitude ranges variates. For this reason, we scale all the data channels to the same amplitude range and
% we normalize each data channel to unit standard deviation so that the total variation becomes equal to the number of channels
% or variables [1].

X=zscore(X');
r = sr*(sum(std(X)));
M = sm*ones(1,size(X,2)); % M: embedding vector
tau=ones(1,size(X,2))*stau; % tau: time lag vector
X=X';


Out_mvMFE=NaN;

for ii=2:Scale
    temp_A=[];
    temp_B=[];
    for iii=1:ii
        Xs = Multi_var(X(:,iii:end),ii);
        [e,A,B] = mvFE(Xs,M,r,n,tau);
        temp_A=[temp_A A];
        temp_B=[temp_B B];
    end
    AA=sum(temp_A);
    BB=sum(temp_B);
    Out_mvMFE(ii)=log(AA/BB);
end
end



function z = Multi_var(Data,to)

%  generate the consecutive coarse-grained time series based on variance
%  Input:   Data: time series;
%           S: the scale factor
% Output:
%           M_Data: the coarse-grained time series at the scale factor S

z=[];
for j=1:size(Data,1)
    m=Data(j,1:floor(size(Data,2)/to)*to);
    y=reshape(m,to,[]); % y: to*[]
    z=[z;var(y)];
end
end