function    [L,dL,ddL]=FitMaxEntLinearCostFunction(lambda,x,means)
%objective function (likelihood) for fitting a discrete maximum entropy (or
%log-linear) model of the form P(x)=1/z exp(lambda'*x) to means "means"

%input: parameters lambda
%matrix of all possible states x of the system
%means: means that we want to match

%output: (normalized) likelihood of data with means "means", and its
%gradient and hessian

[logP,logZ,P]=logPMaxEnt(x,lambda);

L=logZ-means*lambda;

%keyboard
dL= x'*P-means';

if nargout==3
ddL=x'*bsxfun(@times,x,P);
end
