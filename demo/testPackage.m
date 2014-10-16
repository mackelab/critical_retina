%clear all
clc

% Simulation setup
%--------------------------------------------------------------------------
d=25; %simulate 10 dimensional problem
nSamplesTrain = 1000; %[100,1000,10000]; %generate 1000 data-points;
nSamplesEval  = 1000; %[100,1000,10000,100000];
burnIn        = 1000;  %[100,1000,10000,10000];
thinning      = [1,1,1,1];
modelTrue = 'ising_count_l_0';
modelFit  = 'ising_count_l_0';
newLambda = false;
newData   = false;

if newLambda
 h=randn(d,1)-1; %generate random bias terms;
 J=randn(d); J=triu(J,1)/sqrt(d); 
 %J=0*J;
 lambda=hJ2lambda(h,J);
 switch modelTrue
     case 'ising_count_l_0'
         L=randn(d+1,1)/sqrt(d);
     case 'ising_count'
         L=randn(d,1)/sqrt(d);
     case 'ising'
         L = [];
 end
 lambdaTrue = [lambda;L];
end



figure;  
for r = 1:length(nSamplesTrain)
% generate training data
%--------------------------------------------------------------------------
if newData
  disp('Generating training data')
  % Initialize training data-generating MCMC chain with a sample drawn from
  % a nested model (only h = lamdbdaTrue(1:d), i.e. no J, no L)
  EX = exp(lambdaTrue(1:d))./(1+exp(lambdaTrue(1:d))); % works ok for small
  x0 = double(rand(d,1)<EX);                           % entries of L,J

  xTrain = maxEnt_gibbs(nSamplesTrain(r), burnIn(r), thinning(r), ...
                        lambdaTrue, x0, modelTrue, 'default', 'samples');
mfxTrain = xTrain;
%[fxTrain, fDescr] = setup_features_maxent(xTrain', modelFit);
%mfxTrain = mean(fxTrain, 1)';  
 clear fxTrain
 % for comparison: evaluate independent model, i.e. J = 0,L = 0
 EX = mean(xTrain,2);
 idxL = length(lambda)+1:length(lambdaTrue);
 lambdaInd = zeros(size(lambdaTrue)); mfxInd = zeros(size(lambdaTrue));
 lambdaInd(1:d) = log(EX./(1-EX));
 lambdaInd(lambdaInd==Inf) =   1000; % fairly hacky solution in case
 lambdaInd(lambdaInd==-Inf) = -1000; % EX(k) = 0 resp. EX(k) = 1

 mfxInd(1:d) = EX;
 tmp = histc(sum(bsxfun(@lt, EX, rand([d,nSamplesTrain(r)])),1),0:d+1)/nSamplesTrain(r);
 if ~isempty(idxL)
  mfxInd(idxL) =tmp(1:d+1);
 end
 clear EX tmp
end

% train model
%--------------------------------------------------------------------------
fitoptions = struct;
fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
%fitoptions.display='off';
fitoptions.MaxIter=3000;
fitoptions.lambda0 = -lambdaInd;
disp('Fitting maxEnt model')
[lambdaMPF,logZ,logP,fitmeans,output] = fit_maxent_mpf(xTrain',fitoptions);
lambdaMPF = -lambdaMPF;

% validate model
%--------------------------------------------------------------------------
 disp('Generating data from model fit')
% [~, x0] = max(logP); x0 = xTrain(:,x0);
  mfxEval = maxEnt_gibbs_pair(nSamplesEval(r), burnIn(r), thinning(r), ...
                       lambdaMPF, x0, modelFit, 'pair', 'means');
% fxEval = setup_features_maxent(xEval', modelFit);
%mfxEval = mean(fxEval, 1)';
 clear fxEval
% compute true moments
%--------------------------------------------------------------------------
if d < 20
 [features,description,x]=setup_features_maxent(d,modelTrue);
 [logPtrue,logZtrue,Ptrue, means_true]=logPMaxEnt(features,lambdaTrue);
 EX = sum(bsxfun(@times, x', Ptrue'),2);
 description(isnan(description)) = d+1;
 x1 = x; x1(:,end+1) = 1; 
 EXX = sum(bsxfun(@times, (x1(:,description(1,d+1:d*(d+1)/2))...
                        .* x1(:,description(2,d+1:d*(d+1)/2)))',Ptrue'),2);
 EK = zeros(length(L),1);
switch modelTrue
    case 'ising_count'
     for k = 1:length(EK)
      EK(k) = sum((sum(x,2)==(k)) .* Ptrue);
     end
    case 'ising_count_l_0'
     for k = 1:length(EK)
      EK(k) = sum((sum(x,2)==(k-1)) .* Ptrue);
     end
end
mfxTrue = [EX(:);EXX(:);EK(:)];
clear x1 EX EXX EK

end % if d < 20  
 
% depict results
%--------------------------------------------------------------------------
 disp('Visualizing results')
subplot(length(nSamplesTrain),20,(1:15)+(r-1)*20) 
  % subplot for true and recovered moments E[f(X)]
hold on
if d < 20
 plot(mfxTrue,  'r--', 'linewidth', 1.5)   
 plot(mfxTrain, 'ko-', 'linewidth', 1.5) 
 plot(mfxEval,  'bo-', 'linewidth', 1.5) 
 plot(1:d, mfxInd(1:d),   'co-', 'linewidth', 1)
 plot(idxL, mfxInd(idxL),   'co-', 'linewidth', 1)
 axis([0.5, length(mfxEval)+0.5, ...
           min([mfxTrain(:);mfxEval(:)]),max([mfxTrain(:);mfxEval(:)])])
 line([d,d]+            0.5,[min([mfxTrain(:);mfxEval(:)]), ... 
                 max([mfxTrain(:);mfxEval(:)])], 'color', 'k') 
 line([d,d]+(d*(d-1)/2)+0.5,[min([mfxTrain(:);mfxEval(:)]), ... 
                 max([mfxTrain(:);mfxEval(:)])], 'color', 'k') 
 if r == 1
 legend('true', 'training samples', 'model fit samples', ...
        'ind. model (only h)',  'location', 'north')
 end
else
 plot(mfxTrain, 'ko-', 'linewidth', 1.5) % replot these one 
 plot(mfxEval,  'bo-', 'linewidth', 1.5)  
 plot(1:d,mfxInd(1:d),   'co-', 'linewidth', 1)
 plot(idxL, mfxInd(idxL),   'co-', 'linewidth', 1)
 axis([0.5, length(mfxEval)+0.5, min(mfxEval(:)),max(mfxEval(:))])
 line([d,d]+ 0.5,[min(mfxEval(:)),max(mfxEval(:))],'color','k')% delineate 
 line([d,d]+(d*(d-1)/2)+0.5,[min(mfxEval(:)), ...              % parameters
                             max(mfxEval(:))], 'color', 'k')   % for h,J,L
 if r == 1
 legend('training samples', 'model fit samples', 'ind. model (only h)', ...
        'location', 'north')
 end
end
hold off
if r == 1 
 title(['True and recovered distribution moments E[f_i(X)], ', ...
        '#Samples_{train} = ',  num2str(nSamplesTrain(r)), ...
      ', #Samples_{eval} = ',   num2str(nSamplesEval(r)), ...
      ', #Samples_{burnIn} = ', num2str(burnIn(r)),  ...
                          ', ', num2str(thinning(r)), '-thinning'] )    
 
else
 title(['#Samples_{train} = ',  num2str(nSamplesTrain(r)), ...
      ', #Samples_{eval} = ',   num2str(nSamplesEval(r)), ...
      ', #Samples_{burnIn} = ', num2str(burnIn(r)),  ...
                          ', ', num2str(thinning(r)), '-thinning'] )
end        
if r < length(nSamplesTrain)
 set(gca, 'XTick', []);
else
 xlabel('i'); 
end
ylabel('$$E[f_i(X)]$$', 'Interpreter', 'Latex')  
box off; set(gca, 'TickDir', 'out')
legend boxoff

% Add subplot for MSE comparison between different chains and 
% discretely sampled stuff. Always take 'true' moments as reference
lbls = {'full $\lambda$', 'h', 'J', 'L'};

subplot(length(nSamplesTrain),20,(17:20)+(r-1)*20) 
idx = zeros(4,length(mfxTrain));
idx(1,1:length(mfxTrain)) = 1;
idx(2,1:d) = 1; 
idx(3,d+1:d*(d+1)/2) = 1;
idx(4,end-d:end) = 1;
idx = logical(idx);
mse = zeros(4,1); 
mseInd = zeros(4,1);
for i = 1:4
 mse(i) = mean(bsxfun(@minus, mfxTrain(idx(i,:)), mfxEval(idx(i,:))).^2,1);
end
plot(1:4,mse, 'o--', 'linewidth', 2)
if r == length(nSamplesTrain)
 set(gca, 'XTick', 1:4);
 set(gca, 'XTickLabel', lbls);
 rotateticklabel(gca,-45); % double trick: keep stuff readable AND get
                           % LaTeX stuff to be read for 'XTickLabel'
else
 set(gca, 'XTick', [])
end
if r == 1
 title('Mean squared errors')
end
box off
axis([0.5, 4.5, 0.8*min(mse), max(1.1*max(mse),10e-20)]) 
%ylabel('$$\frac{1}{d} \Sigma_{k=1}^d (E[f_k(X)] - \hat{E}[f_k(X)])^2$$',...
%       'Interpreter', 'Latex')

end
% Finalize image by giving it a nice super-title
switch modelFit
    case 'ising'
        modelFit = 'Ising';
    case 'ising_count'
        modelFit = 'Ising + V(K)';
    case 'ising_count_l_0'
        modelFit = 'Ising + V(K)';
end
switch modelTrue
    case 'ising'
        modelTrue = 'Ising';
    case 'ising_count'
        modelTrue = 'Ising + V(K)';
    case 'ising_count_l_0'
        modelTrue = 'Ising + V(K)';
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Model fitting results for fitted model ', modelFit, ...
               ', true model ', modelTrue, ', d = ', ...
                num2str(d)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
