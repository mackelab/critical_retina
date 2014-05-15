function [P,Z]=flat_ising_count_distrib(h,J,n,logit)

if nargin==3
    logit=false;
end

range=[0:n];


hofx=(gammaln(n+1)-gammaln(range+1)-gammaln(n-range+1));

H=hofx+range*h+0.5*range.^2*J;


logZ=logsumexp(H');

H=H-logZ;

if logit
    P=(H);
    Z=logZ;
else
    P=exp(H):
    Z=exp(logZ);
end

