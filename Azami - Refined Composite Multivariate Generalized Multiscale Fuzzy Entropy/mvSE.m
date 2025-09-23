function [Out_mvSE,p1,p2]=mvSE(X,M,r,tau)
%
% This function calculates multivariate sample entropy (mvSE) of a multivariate signal
%
% Inputs:
%
% X: multivariate signal - a matrix of size nvar (the number of channels) x nsamp (the number of sample points for each channel)
% M: embedding vector
% r: scalar threshold
% tau: time lag vector
%
% Outputs:
%
% e: scalar quantity - the mvSE of X
% p1: scalar quantity - the probability that any two composite delay vectors are similar in dimension m
% p2 : scalar quantity - the probability that any two composite delay vectors are similar in dimension m+1
%
% The code is based on the publicly-available code in
% "http://www.commsp.ee.ic.ac.uk/~mandic/research/Complexity_Stuff.htm"
% prepared by Prof. Mandic's group
%
% Ref:
% [1] M. U. Ahmed and D. P. Mandic, "Multivariate multiscale entropy
% analysis", IEEE Signal Processing Letters, vol. 19, no. 2, pp.91-94.2012
% [2] H. Azami and J. Escudero, "Refined Composite Multivariate Generalized Multiscale Fuzzy Entropy:
% A Tool for Complexity Analysis of Multichannel Signals", Physica A, 2016.
%
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
%  10-June-16

mm=max(M);
mtau=max(tau);
nn=mm*mtau;


[nvar,nsamp]=size(X);
N=nsamp-nn;
A=embd(M,tau,X);%all the embedded vectors are created
y=pdist(A,'chebychev');%infinite norm is calculated between all possible pairs
[r1,c1,v1]= find(y<=r);% threshold is implemented
p1=numel(v1)*2/(N*(N-1));%the probability that two templates of length m are closed within the tolerance r
clear  y r1 c1 v1 A;

M=repmat(M,nvar,1);
I=eye(nvar);
M=M+I;

B=[];

% number of match templates of length m+1 closed within the tolerance r where m=sum(M) is calculated afterwards
for h=1:nvar
    Btemp=embd(M(h,:),tau,X);
    B=vertcat(B,Btemp);% all the delay embedded vectors of all the subspaces of dimension m+1 is concatenated into a single matrix
    Btemp=[];
end
z=pdist(B,'chebychev'); %now comparison is done between subspaces
[r2,c2,v2]= find(z<=r);
p2=numel(v2)*2/(nvar*N*(nvar*N-1));
clear  z r2 c2 v2 B;

Out_mvSE=log(p1/p2);
