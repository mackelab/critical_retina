function testPackageCluster(fname,d,nSamplesTrain,nSamplesEval,burnIn, ...
                               nSamplesIterScaling, burnInIterScaling, ...
                               maxIter, maxInnerIter, beta, ...
                               lambda0, lambdaTrue, mfxTrain)
% Simulation setup
%--------------------------------------------------------------------------
folder = '/home/marcel/criticalityIterScaling/results/';
if nargin < 1 || isempty(fname)
 fname = [folder, 'tmp.mat'];
else
 fname = [folder, fname];
end
if nargin < 2 || isempty(d)
d=15; %simulate 15 dimensional problem
end
if nargin < 3 || isempty(nSamplesTrain)
nSamplesTrain = 10000; %[100,1000,10000]; %generate 1000 data-points;
end
if nargin < 4 || isempty(nSamplesEval)
nSamplesEval  = 10000; %[100,1000,10000,100000];
end
if nargin < 5 || isempty(burnIn)
burnIn        = 1000;  %[100,1000,10000,10000];
end
if nargin < 6 || isempty(nSamplesIterScaling)
 fitoptions.nSamples = 10000;    
else 
 fitoptions.nSamples = nSamplesIterScaling;
 clear nSamplesIterScaling;
end
if nargin < 7 || isempty(burnInIterScaling)
 fitoptions.burnIn  = 1000;
else
 fitoptions.burnIn = burnInIterScaling;
 clear burnInIterScaling;
end
if nargin < 8 || isempty(maxIter)
 fitoptions.maxIter = 10*d;
else
 fitoptions.maxIter = maxIter;
 clear maxIter
end
if nargin < 9 || isempty(maxInnerIter)
 fitoptions.maxInnerIter = 1;    
else
 fitoptions.maxInnerIter = maxInnerIter;
 clear maxInnerIter
end
if nargin < 10 || isempty(beta)
 beta = 0.001*ones(d*(d+1)/2 + d+1,1);
end
if nargin < 11
 useLambdaInd = true; 
else 
 fitoptions.lambda0 = lambda0;
 clear lambda0;
end

% The following shall be fixed to always these value
% general sampling/data generation/evaluation related (IT overhead...)
thinning      = ones(size(nSamplesTrain)); % essentially means no thinning
modelTrue = 'ising_count_l_0'; % full model
modelFit  = 'ising_count_l_0'; % full model
trainMethod = 'iterativeScaling';
% iterative scaling-related
fitoptions.regular = 'l1';
fitoptions.machine  = 'cluster';
fitoptions.nRestart = 1;
fitoptions.modelFit = modelFit;

if nargin < 12 || isempty(lambdaTrue)
 h=randn(d,1)-1; %generate random bias terms;
 J=randn(d); J=triu(J,1)/sqrt(d); 
 lambda=hJ2lambda(h,J);
 switch modelTrue
     case 'ising_count_l_0'
         L=randn(d+1,1)/sqrt(d);
     case 'ising_count'
         L=randn(d,1)/sqrt(d);
     case 'ising'
         L = [];
 end
 L = L - L(1); % a rather recently discovered nasty extra degree of freedom
 lambdaTrue = [lambda;L];
end

if nargin < 13 || isempty(mfxTrain)
 newData = true;
end

compTimes = zeros(size(nSamplesTrain));

for r = 1:length(nSamplesTrain)
tic
% generate training data
%--------------------------------------------------------------------------
if newData
  disp('Generating training data')
  % Initialize training data-generating MCMC chain with a sample drawn from
  % a nested model (only h = lamdbdaTrue(1:d), i.e. no J, no L)
  EX = exp(lambdaTrue(1:d))./(1+exp(lambdaTrue(1:d))); % works ok for small
  x0 = double(rand(d,1)<EX);                           % entries of L,J

 [mfxTrain,~,~] = maxEnt_gibbs_pair_C(nSamplesTrain(r), burnIn(r), lambdaTrue, x0, fitoptions.machine);
 %xTrain = maxEnt_gibbs_pair(nSamplesTrain(r), burnIn(r), thinning(r), ...
 %                       lambdaTrue, x0, modelTrue, 'default', 'samples');
 EX = mfxTrain(1:d);
 idxL = length(lambda)+1:length(lambdaTrue);
 lambdaInd = zeros(size(lambdaTrue)); mfxInd = zeros(size(lambdaTrue));
 lambdaInd(1:d) = log(EX./(1-EX));
 lambdaInd(lambdaInd==Inf) =   1000; % fairly hacky solution in case
 lambdaInd(lambdaInd==-Inf) = -1000; % EX(k) = 0 resp. EX(k) = 1

 if useLambdaInd 
   fitoptions.lambda0 = lambdaInd;
 end
 
 mfxInd(1:d) = EX;
 tmp = histc(sum(bsxfun(@lt, EX, rand([d,nSamplesTrain(r)])),1),0:d+1)/nSamplesTrain(r);
 if ~isempty(idxL)
  mfxInd(idxL) =tmp(1:d+1);
 end
 clear EX tmp xTrain fxTrain 
end

% train model
%--------------------------------------------------------------------------
disp('Fitting maxEnt model')
switch trainMethod
  case 'MPF'
    fitoptions = struct;
    fitoptions.optTol=1e-100; 
    fitoptions.progTol=1e-100; 
   %fitoptions.display='off';
    fitoptions.MaxIter=3000;
    fitoptions.lambda0 = -lambdaInd;
    [lambdaHat,logZ,logP,fitmeans,output] = fit_maxent_mpf(xTrain',fitoptions);
    lambdaHat = -lambdaHat;
  case 'iterativeScaling'
    
    %fitoptions.lambda0 = lambdaTrue;
    %fitoptions.lambda0(1:(d*(d+1)/2)) = 0;
    %fitoptions.lambda0(1:d) = lambdaInd(1:d);
    
    %fitoptions.lambda0 = lambdaTrue;
    %fitoptions.lambda0(end-d:end) = 0;
    
    [lambdaHat, fitDiagnostics] = iterScaling(mfxTrain, fitoptions, beta);
end
% validate model
%--------------------------------------------------------------------------
 disp('Generating data from model fit')
 [mfxEval,~,~] = maxEnt_gibbs_pair_C(nSamplesEval(r), burnIn(r), lambdaHat, x0, fitoptions.machine);
%--------------------------------------------------------------------------
% if d < 20
%  [features,description,x]=setup_features_maxent(d,modelTrue);
%  [~,~,Ptrue, ~]=logPMaxEnt(features,lambdaTrue);
%  EX = sum(bsxfun(@times, x', Ptrue'),2);
%  description(isnan(description)) = d+1;
%  x1 = x; 
%  x1(:,end+1) = 1; 
%  EXX = sum(bsxfun(@times, (x1(:,description(1,d+1:d*(d+1)/2))...
%                         .* x1(:,description(2,d+1:d*(d+1)/2)))',Ptrue'),2);
%  EK = zeros(length(L),1);
% switch modelTrue
%     case 'ising_count'
%      for k = 1:length(EK)
%       EK(k) = sum((sum(x,2)==(k)) .* Ptrue);
%      end
%     case 'ising_count_l_0'
%      for k = 1:length(EK)
%       EK(k) = sum((sum(x,2)==(k-1)) .* Ptrue);
%      end
% end
% mfxTrue = [EX(:);EXX(:);EK(:)];
% clear x1 EX EXX EK description features Ptrue 
% end % if d < 20  
compTimes(r) = toc;
end % end for r = 1:length(nSamplesTrain)

if ~exist('mfxTrue', 'var')
    mfxTrue = [];
end

pars.d = d; 
pars.beta = beta; 
pars.nSamplesTrain = nSamplesTrain;
pars.nSamplesEval  = nSamplesEval;
pars.burnIn = burnIn;

save(fname, 'lambdaTrue', 'lambdaHat',          'lambdaInd', ...
            'mfxTrain',   'mfxEval', 'mfxTrue', 'mfxInd', ...
            'fitoptions', 'fitDiagnostics', 'pars', 'compTimes');

end