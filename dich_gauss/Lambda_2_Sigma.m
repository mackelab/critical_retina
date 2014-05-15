function [mu, Rho]=Lambda_2_Rho(gamma, Lambda)
%given latent gaussian, find mean and variance of resulting binary rvs

variances=diag(Lambda);
gamma=gamma./sqrt(variances);
Lambda=cov_2_corr(Lambda);

mu=normcdf(gamma);

Rho=ones(numel(gamma));
for k=1:numel(gamma)
    Rho(k,k)=mu(k).*(1-mu(k));
    for kk=k+1:numel(gamma);
        Rho(k,kk)=bivnor(-gamma(k),-gamma(kk),Lambda(k,kk))-mu(k)*mu(kk);
        Rho(kk,k)=Rho(k,kk);
    end
end

