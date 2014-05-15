function [logP,logZ,P, means]=PMaxEnt(x,lambda,logZ);
%calculate log probabilities, log partition function, probabilites and
%feature means under a model of the form P(x)=exp(-logZ+lambda'*x)

if numel(x)==1
    dimo=x;
    %use standard ordering for x, assuming a second order Ising model:
    %dimo=dimo=-1+sqrt(1+4*numel(lambda));
    [x]=SetupFeaturesMaxEnt(dimo,2);
   % keyboard
end

logQ=full(x*lambda);

if nargin==2
logZ=logsumexp(logQ);
end

logP=logQ-logZ;


if nargout>=3
P=exp(logP);
SumP=sum(P);
P=P/SumP;
end

if nargout==4
    means=P(:)'*x;
end


