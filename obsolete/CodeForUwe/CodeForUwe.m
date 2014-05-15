%A little demo-script for demonstrating the code which attempts to guess
%input-correlations from measured spike-train statistics. 

%For more details and applications, see 
%Dorn, Jessy D., and Dario L. Ringach. "Estimating membrane voltage correlations from extracellular spike trains." Journal of neurophysiology 89.4 (2003): 2271-2278.
%Macke, Jakob H., et al. "Generating spike trains with specified correlation coefficients." Neural Computation 21.2 (2009): 397-423.
%Macke, Jakob H., Manfred Opper, and Matthias Bethge. "Common input explains higher-order correlations and entropy in a simple model of neural population activity." Physical review letters 106.20 (2011):

%We will assume that your data is from binary spikes, i.e. that the first
%neuron spikes with probability p_1 and the second with probability p_2,
%and that the (Pearson) correlation coefficient between the two outputs is
%given by rho. The code will calculate the input-correlation between the
%two neurons assuming that each neuron is a simple threshold device, and
%that their inputs are Gaussian. In principle, it is easy to do the
%corresponding calculations for more complicated models (and, in
%particular, ones which are matched to the whole-cell data), but this would
%require some more work.

%try this out on some fake-data, 5 pairs of neurons:
p_1=[0.1, 0.02,  .1,  .1, .02, .1,   .02];
p_2=[0.2, .01, .02, .03, .1, .04, .02];
rho=[.4,  .03, .05,  .3, .4, .1,  .1];

[lambda]=FindDGcorrelation(p_1, p_2, rho)



plot(rho,lambda,'.')