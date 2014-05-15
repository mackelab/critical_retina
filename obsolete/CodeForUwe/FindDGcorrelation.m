function [lambda]=FindDGcorrelation(p_1, p_2, rho)
%calculate the 'estimated input correlation' lambda assuming a threshold
%neuron with Gaussian inputs (a.k.a a "Dichotomized Gaussian"). 

%Inputs: 
%p_1, an n by 1 vector of the firing probaility of the first neuron
%in each pair
%p_2, an n by 1 vector of the firing probability of the second neuron
%rho%, an n by 1 vector of the Pearson correlation coefficient of the two
%neurons
%output: 
%lambda, an n by 1 vector of the estimated input correlations
%
%JHM 01/2013

%For more details and applications, see 
%Dorn, Jessy D., and Dario L. Ringach. "Estimating membrane voltage correlations from extracellular spike trains." Journal of neurophysiology 89.4 (2003): 2271-2278.
%Macke, Jakob H., et al. "Generating spike trains with specified correlation coefficients." Neural Computation 21.2 (2009): 397-423.
%Macke, Jakob H., Manfred Opper, and Matthias Bethge. "Common input explains higher-order correlations and entropy in a simple model of neural population activity." Physical review letters 106.20 (2011):


n=numel(p_1);

for k=1:n
    %mean vector:
    mu=[p_1(k); p_2(k)];
    %construct covariance from correlations:
    varo=mu.*(1-mu);
    c=rho(k)*sqrt(prod(varo));
    Sigma=[varo(1), c; c, varo(2)];
    try
    [gammaloc,lambdaloc]=findLatentGaussian(mu,Sigma);
    lambda(k)=lambdaloc(2);
    catch
        lambda(k)=nan;
        warning('An error occured while computing input correlation-- result replaced by nan')
    end
end


