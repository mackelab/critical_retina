function [h,J,logZ,logP, patterns]=FitIsingModel(mu,Cov)
%function [h,J,logZ,logP, patterns]=FitIsingModel(mu,Cov)
%
%fits parameters of an Ising model (i.e. second order binary maximum entropy model)
% to binary data with mean mu and covariance Cov. Assumes that data is
% represented at {0,1} (not {-1,1}). In other words, finds a vector h and
% matrix J such that the distribution P(x)=1/z exp(h'x+ 0.5*x' J x) has mean mu
% and covarianve Cov.
%
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

%uses minFunc by Mark Schmidt


%find dimensionality of input space
d=numel(mu);


%need to convert mu and Cov to the binary feature expectations that I
%usually work with (i.e. P(x)=1/z exp (lambda' *features(x)))
%make all patterns, and corresponding feature-representations of all
%2-tupels on d binary patterns:
[features,description,patterns]=SetupFeaturesMaxEnt(d,2);


%calculate mean of second order features, i.e. second order
%correlation matrix
triup=triu(ones(size(Cov)),1)==1;
E2=Cov+mu(:)*mu(:)';

pairmeans=E2(triup);

%get overall feature expecations by concatenating means and upper triangle
%of correlation matrix:
means=[mu(:);pairmeans]';

%use general purpose function "FitMaxEntLinear" to learn parameters:
fitoptions.TolX=1e-20;
fitoptions.TolFun=1e-20;
fitoptions.display='off';
[lambda,logZ, logP, junk,junk2]=FitMaxEntLinear(features,means, fitoptions);

%now, extract h and J from the weights lambda:
[h,J]=hJ2lambda(lamdba);
