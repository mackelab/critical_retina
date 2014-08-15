function [h,J,logZ,logP, patterns,l]=fit_ising_model(mu,Cov,PK)
%function [h,J,logZ,logP, patterns,l]=fit_ising_model(mu,Cov)
%
%fits parameters of an Ising model (i.e. second order binary maximum entropy model)
% to binary data with mean mu and covariance Cov. Assumes that data is
% represented at {0,1} (not {-1,1}). In other words, finds a vector h and
% matrix J such that the distribution P(x)=1/z exp(h'x+ 0.5*x' J x) has mean mu
% and covarianve Cov.
%
%If a third argument, PK is used, then it incoproprates the additional
% constraints that the activity distribution PK= P(sum(x)=k) is matched, as
% in Tkacik et al, arxiv 2014. PK should be a vector of the same size as
% mu, with PK(k)= P(sum(x)=k), then the model being fit is
% P(x)=1/z exp(h'x+ 0.5*x' J x+ l_k sum(x))

%inputs:
%mu: vector of mean activities
%Cov: covariance of activities
%
%outputs:
%h:  vector of bias terms h
%J: matrix of interaction terms J
%logZ: log partition function, i.e. log of normalizer of distribution
%logP: for each possible binary pattern, its log-probability
%patterns: a vector all binary patterns with as many elements as mu

%l: (only if PK is used) vector of weights for spike-count distribution


%uses minFunc by Mark Schmidt


%find dimensionality of input space
d=numel(mu);
if nargin==3
    count_constraints=true;
    if numel(PK)>d
        error('PK should only contain probabilities for counts>0');
    end
else
    count_constraints=false;
end

%need to convert mu and Cov to the binary feature expectations that I
%usually work with (i.e. P(x)=1/z exp (lambda' *features(x)))
%make all patterns, and corresponding feature-representations of all
%2-tupels on d binary patterns:
if ~count_constraints
    [features,description,patterns]=setup_features_maxent(d,2);
    means=meancov_2_features(mu,Cov);
    penalties=ones(size(means'))*1e-9;
else
    [features,description,patterns]=setup_features_maxent(d,'ising_count');
    means=meancov_2_features(mu,Cov);
    means=[means, PK];
    penalties=ones(size(means'))*1e-9;
    penalties(end-d+1:end)=1e-6;
end

%get overall feature expecations by concatenating means and upper triangle
%of correlation matrix:
%means=[mu(:);pairmeans]';

%use general purpose function "FitMaxEntLinear" to learn parameters:
fitoptions.optTol=1e-20;
fitoptions.progTol=1e-20;
fitoptions.display='off';
[lambda,logZ, logP, junk,junk2]=fit_maxent_linear(features,means, fitoptions,0,penalties);

%now, extract h and J from the weights lambda:
if ~count_constraints
    [h,J]=hJ2lambda(lambda);
else
    [h,J]=hJ2lambda(lambda(1:end-d));
    l=lambda(end-d+1:end);
end
