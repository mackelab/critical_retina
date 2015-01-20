function getEfxEval(idxRun, nSamples, burnIn)
% computes one fairly large last sample from the last known parameter
% vector. 

if nargin < 3
    burnIn = 512000;
end
if nargin < 2
    nSamples = 512000;
end
for i = idxRun
 
 load(['/home/marcel/criticalityIterScaling/s', num2str(i),'_res_small.mat'])
 addpath(genpath('/home/marcel/criticalityIterScaling/code/'))
 x0 = 100; % urges the sampler to randomly initialize the first draw
           % from the independent model
 disp(['will sample from last lambda of run ', num2str(i)])
 disp(['nSamples = ', num2str(nSamples)])
 disp(['burnIn = ', num2str(burnIn)])
 disp('Last chance to abort, as sampler is C++ and cannot be ctrl+c-killed')
 disp('Press key to continue')
 pause;
 disp('Ok, starting to do sampling!')           
[out.EfxHat,~,~] = maxEnt_gibbs_pair_C(nSamples, burnIn, out.lambdaHat, x0, 'cluster');

 disp('done sampling, storing in s*_res_small.mat')
 save(['/home/marcel/criticalityIterScaling/s', num2str(i),'_res_small.mat'], 'out');

 disp('storing in EfxHat.mat')
 load('/home/marcel/criticalityIterScaling/EfxHat.mat')
 EfxHat(:,i) = out.EfxHat; 
 save('/home/marcel/criticalityIterScaling/EfxHat.mat', 'EfxHat')
 
end
end