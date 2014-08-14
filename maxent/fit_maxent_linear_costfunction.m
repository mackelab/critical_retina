function    [L,dL,ddL]=fit_maxent_linear_costfunction(lambda,x,means,weights)
%objective function (likelihood) for fitting a discrete maximum entropy (or
%log-linear) model of the form P(x)=1/z exp(lambda'*x) to means "means"

%input: parameters lambda
%matrix of all possible states x of the system
%means: means that we want to match
%weights: set of fixed weights to be added

%output: (normalized) likelihood of data with means "means", and its
%gradient and hessian

[logP,logZ,P]=logPMaxEnt(x,lambda,[],weights);

L=logZ-means*lambda-mean(weights);

%keyboard
dL= x'*P-means';

if nargout==3
ddL=x'*bsxfun(@times,x,P);
end
