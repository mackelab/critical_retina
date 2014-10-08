function [out] = retSim(x, W, Cn, mode)
% input:
%  -  x: d-by-N vector, gives a sequence of N input images
%  -  W: n-by-d matrix of linear filters for each neuron i = 1,...,n
%  - Cn: n-by-n covariance matrix for noise correlations

[n,d] = size(W);
[~,N] = size(x);

y = W * x; % linear filtering for generating RGC input
e = mvnrnd(zeros(1,n), Cn, N)'; % generate noise
y = y + e; % add noise on level of RGC input

y = 1 ./ (1 + exp(-y)); % compute RGC output firing probability

out.spikes = (rand([n,N])<y);

end