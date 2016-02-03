function [cN, pcountT] = computeSpecificHeatTraces(pcount, Ts, ifApproxByBetaBinomial, alphabeta)
% function version to compute specific heat traces for flat models from 
% (empirical) P(K). Has the option to approximate P(K) with a beta-binomial
% distribution, which allows predicting specific heat at T=1 for large n.

if nargin < 2    
    Ts = (0.8 : 0.0125 : 2)'; % range of temperatures to look at 
end
n = length(pcount)-1;                
lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';

if nargin < 3
    ifApproxByBetaBinomial = true; 
end

if ifApproxByBetaBinomial
    if nargin < 4
        mu1 =  (0:n) * pcount;
        mu2 = (0:n).^2 * pcount;
        Z = ( n * (mu2/mu1 - mu1 -1)) + mu1;
        a = (n * mu1 - mu2) / Z;
        b = (n - mu1) * (n - mu2/mu1) / Z;
    else
        a = alphabeta(1);
        b  = alphabeta(2);
    end
    logpcount = lognchoosek + betaln(a + (0:n), n+b-(0:n))'-betaln(a,b);
    pcount = exp(logpcount);
end

lambdaT1 = zeros(n+1,1);                     % get lambda from count distr.
lambdaT1(2:n+1) = fit_flat_maxent_model(pcount); % fix lambda_0 = 0 

pcountT = zeros(n+1,length(Ts));
varE    = zeros(length(Ts),1);
cN      = zeros(length(Ts),1);
for i = 1:length(Ts)
    T = Ts(i);
    lambdaT = lambdaT1/T;
    
    logpcountT = -Inf*ones(n+1,1);
    logpcountT(1:n+1) = lambdaT(1:n+1) + lognchoosek;
    pcountT(:,i) = exp(logpcountT);
    pcountT(:,i) = pcountT(:,i)/sum(pcountT(:,i)); % getting Z is easy here
    tmp = lambdaT;
    tmp(pcountT(:,i)==0) = 0;  % ensure 0 * Inf = 0
    varE(i) =  (pcountT(:,i)' * tmp.^2) - (pcountT(:,i)' * tmp)^2;
    cN(i) = (1/n) * varE(i); % specific heat, up to Boltzmann constant k_B
end
end



