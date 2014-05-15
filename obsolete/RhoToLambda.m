function [gamma, lambda]=RhoToLambda(mu, rho, usecov, acc)

v(1)=mu(1)*(1-mu(1));
v(2)=mu(end)*(1-mu(end));
if nargin==2 | usecov==0;
    c=rho*sqrt(prod(v));
else
    c=rho;
end
if nargin<=3
    acc=1e-10;
end

covo=[v(1), c; c, v(2)];
[g,L] = findLatentGaussian([mu(1), mu(end)],covo,acc);
gamma=g;
lambda=L(2);
%keyboard

if numel(mu)==1
    gamma=gamma(1);
end
