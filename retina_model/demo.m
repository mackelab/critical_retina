%% ------------------------------------------------------------------------
% (1) Set up the simulation, fitting & evaluation processing chain
%%-------------------------------------------------------------------------
clear all
close all
clc

% Settings concerning the overall retina simulation
%--------------------------------------------------------------------------
d = [20, 20]; % dimension of visual field, in pixels
n = 20;       % number of RGCs to be simulated
pRet.Ce=eye(n); % covariance matrix for Gaussian-induced noise correlations
pRet.magnitude = 1;      % parameters governing the nonlinearity mapping
pRet.gain      = 1;      % from linear filter responses to RGC output
pRet.offset    = -2.944;  % spiking probability
mode.tiling = 'random';           % arrangement of RGC centres 
mode.RF = 'center-surround DoG';  % layout of RGC receptive fields

% Settings concerning input images to be fed into the model retina
%--------------------------------------------------------------------------
N = 10000; % number of image frames to be generated
alpha = 2; % 0 'should' be white, 1 pink and 2 brown noise 
xMag  = 1000; % magnitude of Gaussian noise to generate images from

% Settings concerning the maxEnt model fitting and evaluation
%--------------------------------------------------------------------------
modelFit = 'ising_count_l_0'; % maxEnt model to be fit to data
nSamplesEval = 5000; % number of Gibbs samples to extimate means E[f(X)]
burnIn       =  1000; % number of first Gibbs samples to be discarded
thinning     =     1; % distance in sequence between Gibbs samples to be
                      % stored (integers >1 thin out the MCMC chain)
MLhackData   = true;

% Parameters concerning the specific layout of the RGC receptive fields
%--------------------------------------------------------------------------
pars.thres = 0.00001;% threshold below which the RFs will be truncated to 0
Sigma = { [15,   0 ;      % template for central part of DoG filters
            0,  15 ], ... % (peak of Mexican heat)
          [20,   0 ;      % template for outer part of DoG filters
            0,  20 ]};    % (surround of Mexican hat)
SigmaDFs = 100000*[1,1];  % degrees of freeom for DoG component covariances
ON_OFF = 2*randi(2, [n,1]) - 3; % per cell whether it's an ON or an OFF RGC
hight = [0.9, 1]; % template for hights of central and outer DoG components
hightSTDs = [0.01, 0.01]; % standard deviates for hights of DoG components 
idxON  = find(ON_OFF > 0); % quick lists for finding
idxOFF = find(ON_OFF < 0); % ON and OFF cells
pars.hight = zeros(n,2); % hights of DoG component Gaussian bumps
pars.hight(:,1) = abs(normrnd(hight(1),hightSTDs(1)^2,[n,1])) .*ON_OFF;
pars.hight(:,2) = abs(normrnd(hight(2),hightSTDs(2)^2,[n,1])) .*ON_OFF;
for i = 1:n 
  pars.Sigma{i,1} = wishrnd(Sigma{1}, SigmaDFs(1))/SigmaDFs(1);
  pars.Sigma{i,2} = wishrnd(Sigma{2}, SigmaDFs(2))/SigmaDFs(2);
end

%%------------------------------------------------------------------------
% (2) Generate input images and do the retina simulation
%%-------------------------------------------------------------------------
disp('Generating RGC filters')
[W, RGCcen] = genFilters(d,n,mode,pars); 
W=sparse(W);
%--------------------------------------------------------------------------
disp('Visualizing RGC filters')
figure(1); % 3D images showing individual filters with different colors
subplot(1,2,1)
for i = 1:length(idxON); % ON cells 
 xy = reshape(W(idxON(i),:), d(1), d(2)); 
 xy(xy==0) = NaN;                  % don't plot these pixels, just clutters
 surf(xy, ones(d)*i/length(idxON)), hold on, 
end;
hold off, title('ON RGCs')
subplot(1,2,2)
for i = 1:length(idxOFF); % OFF cells
 xy = reshape(W(idxOFF(i),:), d(1), d(2)); 
 xy(xy==0) = NaN;                  % don't plot these pixels, just clutters
 surf(xy, ones(d)*i/length(idxOFF)), hold on, 
end;
hold off, title('OFF RGCs')
figure(2); % 2D image showing sum of filters, i.e. overall retina imprint
subplot(1,2,1), xy = reshape(sum(W(idxON,:),1), d(1), d(2));
imagesc(xy), title('ON RGCs')
subplot(1,2,2), xy = reshape(sum(W(idxOFF,:),1), d(1), d(2));
imagesc(xy), title('OFF RGCs')
clear xy
%--------------------------------------------------------------------------
disp('Generating input images')
x = xMag * spatialPattern([d(1),d(2),N], -alpha);
%--------------------------------------------------------------------------
disp('Simulating neural activity')
out = retSim(x,W,pRet);
out.spkCorrs = full(corr(out.spikes'));
out.spkCov   = full(cov(out.spikes'));
out.RFoverlap = full(W*W');
%--------------------------------------------------------------------------
disp('Visualizing neural activity')
figure(3); 
subplot(2,3,1:3), imagesc(out.spikes), title('spike trains')
xlabel('time t'), ylabel('neuron n')
set(gca, 'tickdir', 'out'), box off;
subplot(234), imagesc(out.spkCorrs-diag(diag(out.spkCorrs)))
title('output spike correlations')
xlabel('neuron n'), ylabel('neuron n')
set(gca, 'tickdir', 'out'), box off;
subplot(235), imagesc(out.RFoverlap-diag(diag(out.RFoverlap)))
title('receptive field overlaps')
xlabel('neuron n'), ylabel('neuron n')
set(gca, 'tickdir', 'out'), box off;
idxM = (logical(triu(ones(n,n),1))); OnOff = logical((ON_OFF*ON_OFF'+1)/2);
subplot(236), 
plot(out.spkCorrs(idxM&OnOff), out.RFoverlap(idxM&OnOff),'r.');
hold on
plot(out.spkCorrs(idxM&~OnOff), out.RFoverlap(idxM&~OnOff),'g.');
title('spike correlations = f(RF overlap)?')
xlabel('receptive field overlap')
ylabel('output spike correlations')
legend('ON-ON or OFF-OFF', 'ON-OFF or OFF-ON', 'location', 'southeast')
set(gca, 'tickdir', 'out'), box off;

%% ------------------------------------------------------------------------
% (3) Fit the statistical model to data and evaluate the quality of fit
%%-------------------------------------------------------------------------
disp('Fitting maxEnt model')
if MLhackData % current method to avoid infinities in the ML fitting
 out.spikes(:,end+1:end+n) = triu(ones(n,n));
 out.spikes(:,end+1)       = zeros(n,1);
end
% for comparison: evaluate independent model, i.e. J = 0,L = 0
EX = full(mean(out.spikes,2));
idxL = n*(n+1)/2 + (1:n+1);
lambdaInd = zeros(n*(n+3)/2+1,1); 
mfxInd = zeros(size(lambdaInd));
lambdaInd(1:n) = log(EX./(1-EX));
lambdaInd(lambdaInd==Inf) =   1000; % fairly hacky solution in case
lambdaInd(lambdaInd==-Inf) = -1000; % EX(k) = 0 resp. EX(k) = 1
mfxInd(1:n) = EX;
tmp = histc(sum(bsxfun(@lt, EX, rand([n,N])),1),0:n+1)/N;
if ~isempty(idxL)
  mfxInd(idxL) =tmp(1:n+1);
end
clear EX tmp
fitoptions = struct;
fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
fitoptions.display='off';
fitoptions.MaxIter=3000;
fitoptions.lambda0 = -lambdaInd;
[lambdaMPF,~,~,~,~] = fit_maxent_mpf(out.spikes',fitoptions);
lambdaMPF = -lambdaMPF; % conventions...
%--------------------------------------------------------------------------
disp('Generating data from fitted model for evaluation')
 mfxEval = maxEnt_gibbs_pair(nSamplesEval, burnIn, thinning, ...
                       lambdaMPF, n, modelFit, 'rb', 'means');
[fxTrain, fDescr] = setup_features_maxent(out.spikes', modelFit);
 mfxTrain = mean(fxTrain, 1)';  
%--------------------------------------------------------------------------                   
disp('Visualizing fitting results')
figure(4);
subplot(1,20,(1:15))  % subplot for true and recovered moments E[f(X)]
hold on
plot(mfxTrain, 'ko-', 'linewidth', 1.5) 
plot(mfxEval,  'bo-', 'linewidth', 1.5) 
plot(1:n, mfxInd(1:n),   'co-', 'linewidth', 1)
plot(idxL, mfxInd(idxL),   'co-', 'linewidth', 1)
axis([0.5, length(mfxEval)+0.5, ...
          min([mfxTrain(:);mfxEval(:)]),max([mfxTrain(:);mfxEval(:)])])
line([n,n]+            0.5,[min([mfxTrain(:);mfxEval(:)]), ... 
                max([mfxTrain(:);mfxEval(:)])], 'color', 'k') 
line([n,n]+(n*(n-1)/2)+0.5,[min([mfxTrain(:);mfxEval(:)]), ... 
                max([mfxTrain(:);mfxEval(:)])], 'color', 'k') 
legend('training samples', 'model fit samples', ...
       'ind. model (only h)',  'location', 'north')
legend boxoff
title(['True and recovered distribution moments E[f_i(X)], ', ...
        '#Samples_{train} = ',  num2str(N), ...
      ', #Samples_{eval} = ',   num2str(nSamplesEval), ...
      ', #Samples_{burnIn} = ', num2str(burnIn),  ...
                          ', ', num2str(thinning), '-thinning'] )     
xlabel('i'); ylabel('$$E[f_i(X)]$$', 'Interpreter', 'Latex')  
box off; set(gca, 'XTick', []); set(gca, 'TickDir', 'out'), hold off
% Add subplot for MSE comparison between different chains and 
% discretely sampled stuff. Always take 'true' moments as reference
lbls = {'full $\lambda$', 'h', 'J', 'L'};
subplot(1,20,(17:20)) 
idx = zeros(4,length(mfxTrain));
idx(1,1:length(mfxTrain)) = 1;     idx(2,1:n) = 1; 
idx(3,n+1:n*(n+1)/2) = 1;          idx(4,end-n:end) = 1;
idx = logical(idx);
mse = zeros(4,1);  mseInd = zeros(4,1);
for i = 1:4
 mse(i) = mean(bsxfun(@minus, mfxTrain(idx(i,:)), mfxEval(idx(i,:))).^2,1);
end
plot(1:4,mse, 'o--', 'linewidth', 2)
set(gca, 'XTick', 1:4);
set(gca, 'XTickLabel', lbls);
rotateticklabel(gca,-45); % double trick: keep stuff readable AND get
                          % LaTeX stuff to be read for 'XTickLabel'
set(gca, 'XTick', []), title('Mean squared errors'), box off, hold off
axis([0.5, 4.5, 0.8*min(mse), max(1.1*max(mse),10e-20)]) 
% Finalize image by giving it a nice super-title
switch modelFit
    case 'ising'
        modelFit = 'Ising';
    case 'ising_count'
        modelFit = 'Ising + V(K)';
    case 'ising_count_l_0'
        modelFit = 'Ising + V(K)';
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Model fitting results for fitted model ', modelFit, ...
               'n = ', num2str(n)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

figure(5);
pairs=nchoosek(1:n,2); % needed to quickly compute the features of data x
Ceval = zeros(n,n);
for i = 1:size(pairs,1)
   k = pairs(i,1); l = pairs(i,2);
   Ceval(k,l) = mfxEval(i+n)-mfxEval(k)*mfxEval(l);
end
for i = 1:n
   Ceval(i,i) = mfxEval(i)*(1-mfxEval(i));
end
clear k l i
subplot(221), imagesc(Ceval+Ceval'-diag(diag(Ceval)))
axis square, box off, set(gca, 'tickDir', 'out')
xlabel('neuron n'),ylabel('neuron n'), title('Recovered covariance matrix')
subplot(222), imagesc(out.spkCov)
axis square, box off, set(gca, 'tickDir', 'out')
xlabel('neuron n'),ylabel('neuron n'), title('Data covariance matrix')
subplot(223), 
plot(diag(out.spkCov), diag(Ceval), 'r.')
hold on
plot(out.spkCov(idxM), Ceval(idxM), '.')
m = min([Ceval(:);out.spkCov(:)]); M = max([Ceval(:);out.spkCov(:)]);
line([m,M],[m,M], 'color', 'k')
plot(diag(out.spkCov), diag(Ceval), 'r.')
plot(out.spkCov(idxM), Ceval(idxM), '.')
hold off
axis([m,M,m,M])
axis square, box off, set(gca, 'tickDir', 'out')
ylabel('Recovered covariances'),xlabel('Data covariances'), 
title('Recovered vs. data covariance matrix entries')
legend('Var[x_k]', 'Cov[x_k, x_l]','location','southeast')
subplot(224)
plot(0:n,mfxTrain(end-n:end), 'ko-', 'linewidth', 1.5)
hold on
plot(0:n,mfxEval(end-n:end), 'bo-', 'linewidth', 1.5)
xlabel('population activity count K')
ylabel('Pr(K)')
title('Recovered and true population activity counts')
legend('Data activity', 'Recovered activity')
box off, set(gca,'tickDir', 'out'), 
M = max([mfxEval(end-n:end);mfxTrain(end-n:end)]);
axis([-0.5,n+0.5,0,1.05*M])

figure(6); 
plot(lambdaMPF, 'bo-', 'linewidth', 1.5)
line([n,n]+            0.5,[-6,6],'color', 'k')
line([n,n]+ n*(n-1)/2 + 0.5,[-6,6],'color', 'k')
ylabel('\lambda_i'), xlabel('i'), box off
set(gca, 'tickDir', 'out'), title('Model parameter fit')