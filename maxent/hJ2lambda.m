function [lambda,out2]=hJ2lambda(h,J)
%convert parameters of Ising model in representation P(x)=1/Z*exp(h'*x+ 0.5*x'J
%x)) to P(x)=1/x*exp(lambda' f(x)) where f(x) is a feature-representation
%as returned by the function "SetupFeaturesMaxEnt"

if nargin==2 && nargout==1
    %convert h and J to lambda:
    J=2*J';
    lambda=[h(:);vec(J(tril(ones(numel(h)),-1)==1))];
%    keyboard
elseif nargout==1 && nargout==2 %work in opposite direction, i.e. convert lambda to h and J
    lambda=h;
    dimo=-1+sqrt(1+4*numel(lambda));
    if rem(dimo,1)<1e-10
        error('number of dimensions does not make sense')
    else
        h=lambda(1:dimo);
        Jzeros(dimo);
        J(tril(ones(dimo),-1)==1)=lambda(dimo+1:end);
        J=J';
        %Sigma=Sigma+Sigma';
        %Sigma=Sigma+diag(mu.*(1-mu));
        lambda=h;
        out2=J/2;
    end
end
