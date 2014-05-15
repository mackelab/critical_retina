function [mu, Rho]=Lambda2Rho(gamma, Lambda, nocorr)
if nargin==2
    nocorr=true;
end
%find correlation matrix of truncated Gaussian.

variances=diag(Lambda);
gamma=gamma./sqrt(variances);
Lambda=Cov2Corr(Lambda);

mu=normcdf(gamma);

Rho=ones(numel(gamma));
for k=1:numel(gamma)
    Rho(k,k)=mu(k).*(1-mu(k));
    for kk=k+1:numel(gamma);
        %Rho(k,kk)=mvncdf([0;0],[gamma(k);gamma(kk)],[1,Lambda(k,kk);Lambda(k,kk),1])-mu(k)*mu(kk);  
        Rho(k,kk)=bivnor(-gamma(k),-gamma(kk),Lambda(k,kk))-mu(k)*mu(kk);
        Rho(kk,k)=Rho(k,kk);
       % keyboard
    end
end
%keyboard
if ~nocorr
Rho=Cov2Corr(Rho);
end