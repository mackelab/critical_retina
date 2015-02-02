function getEfxEval(idxRun, nSamples, burnIn, timePoint)
% computes one fairly large last sample from the last known parameter
% vector. 

if nargin < 4
    timePoint = 1; % get results after 1 day
end
if nargin < 3
    burnIn = 512000;
end
if nargin < 2
    nSamples = 512000;
end

load('/home/marcel/criticalityIterScaling/timeStampsClusterRuns.mat')
times  = cell(length(idxRun),1);
fnames = cell(length(idxRun),1);
idxTime = zeros(length(idxRun),1); % single chosen point for overview plots

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
 
  tmp = dates(idxSess{i});
  tmp = tmp(2:end, :);   % remove sess*_***.mat, i.e. base file.
  times{i} = datenum(tmp(1:10:end, :));% thin out to match s*_res_small.mat
  tmp = names(idxSess{i});
  tmp = tmp(2:end, :);   % remove sess*_***.mat, i.e. base file.
  fnames{i} = tmp(1:10:end, :);        % thin out to match s*_res_small.mat
  times{i} = times{i} - times{i}(1,:);
  [~, idxTime(i)] = min(abs(times{i}-timePoint));     
  disp(['Extracted time point at ',num2str(timePoint), ...
       ' days it at iteration #', num2str(idxTime(i)),'/',num2str(size(out.lambda,2))])
  lambdaHat = out.lambda(:,idxTime(i));
[out.EfxHat,~,~] = maxEnt_gibbs_pair_C(nSamples, burnIn, lambdaHat, x0, 'cluster');

 %disp('done sampling, storing in s*_res_small.mat')
 %save(['/home/marcel/criticalityIterScaling/s', num2str(i),'_res_small.mat'], 'out');

 disp('storing in EfxHat.mat')
 load(['/home/marcel/criticalityIterScaling/EfxHat',num2str(timePoint),'.mat'])
 EfxHat(:,i) = out.EfxHat; 
 save(['/home/marcel/criticalityIterScaling/EfxHat',num2str(timePoint),'.mat'])
 
end
end