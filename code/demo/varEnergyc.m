
function varEnergyc(n, idxRepet, nSamplesBatch, idxT, Ts)
% This serves to systematically test the length of MCMC samples needed to
% achieve stable estimates of the varience of energy.
% might want to discard the very first estimate Efy(:,1), Moe(:,1), as they
% might have been obtained without any burnin.

if nargin < 2
    idxRepet = 1;
end
if nargin < 3
    nSamplesBatch = 10000;
end
if nargin < 4
    idxT = 1; 
end
if nargin < 5
    %Ts = [0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5];
    Ts = [   0.8000, 0.8563, 0.9000, 0.9363, 0.9675, 0.9938, 1.0175, 1.0388, 1.0575, 1.0775, 1.0988, 1.1200, ...
                1.1413    , 1.1650    , 1.1900    , 1.2150    , 1.2425    , 1.2713    , 1.3013    , 1.3338    , 1.3687    , 1.4063, ...
                1.4463    , 1.4900    , 1.5375 , 1.5900    , 1.6500, 1.7175    , 1.7950 , 1.8875, 2.0000];
end

T = Ts(idxT); 

burnIn = 0;

ticTime = now;
 load('/home/marcel/criticalityIterScaling/results/final/lambdaHatCB2.mat') % got to transfer a couple of lambdas to the cluster,
  % gives 1x10 cell 'lambdaHatCB2'               % ideally those fit to the full model (h,J,V) so far
if idxRepet <= size(lambdaHatCB2{n},2)   
 lambda = lambdaHatCB2{n}(:,idxRepet);
else
 warning('\lambda for desired idxRepet not found - using first one found')
 tmp = lambdaHatCB2{n}(:,1);
 for j = 2:size(lambdaHatCB2{n},2)
  if max(abs(lambdaHatCB2{n}(:,j))) > 0
    tmp = lambdaHatCB2{n}(:,j);
  end
 end
 lambda = tmp; clear tmp
end

n = 10*n;  % move from file index to population size 

nSamplesMax = 1000000000; 
N = floor(nSamplesMax/nSamplesBatch);    %number of update batches                                        
MoE = zeros(2,N); % moments of energy 
VarE = zeros(1,N); 

EX = exp(lambda(1:n))./(1+exp(lambda(1:n))); % a maxEnt model with only h,
x0 = double(rand(n,1)<EX);                   % i.e. no parameters J, L

lambdaT = lambda/T; 
for i = 1:N
   tocTime = now;
   if (tocTime - ticTime) > 1/6 * (n/100)^2;
      break
   end
   disp([num2str(i),'/',num2str(N)])
   [~, MoE(:,i), x0] = maxEnt_gibbs_pair_C(nSamplesBatch, burnIn, lambdaT, x0, 'cluster');
   if mod(i,20)==0
    tmp = num2str(T); tmp(tmp=='.')='_'; 
    save(['/home/marcel/criticalityIterScaling/results/energy/VarEnc', num2str(n), ...
          'idxRepet', num2str(idxRepet), 'T', tmp, '.mat'], ...
          'MoE', 'x0', 'lambdaT', 'nSamplesBatch', 'burnIn')
   end
end
   MoE = MoE(:,1:i);
   tmp = num2str(T); tmp(tmp=='.')='_'; 
    save(['/home/marcel/criticalityIterScaling/results/energy/VarEnc', num2str(n), ...
          'idxRepet', num2str(idxRepet), 'T', tmp, '.mat'], ...
          'MoE', 'x0', 'lambdaT', 'Ts', 'T', 'nSamplesBatch', 'burnIn')
end
