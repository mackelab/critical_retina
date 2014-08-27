clear all
clc

% Simulation setup
%--------------------------------------------------------------------------
d=5; %simulate 10 dimensional problem
nSamplesTrain = 100000; %generate 1000 data-points;
nSamplesEval  = 100000;
burnIn        = 1000;
modelTrue = 'ising_count_l_0';
modelFit  = 'ising_count_l_0';
newLambda = true;
newData   = true;

if newLambda
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
 lambdaTrue = [lambda;L];
end




% generate training data
%--------------------------------------------------------------------------
if newData
  disp('Generating training data')
  x0 = randi(2, [d,1])-1;
  xTrain = maxEnt_gibbs(nSamplesTrain, burnIn, lambdaTrue, x0, modelTrue);
[fxTrain, fDescr] = setup_features_maxent(xTrain', modelTrue);
 fxTrain = fxTrain';
end

% train model
%--------------------------------------------------------------------------
fitoptions = struct;
fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
%fitoptions.display='off';
fitoptions.MaxIter=3000;
disp('Fitting maxEnt model')
[lambdaMPF,logZ,logP,fitmeans,output] = fit_maxent_mpf(xTrain',fitoptions);





% validate model
%--------------------------------------------------------------------------
 disp('Generating data from model fit')
 [~, x0] = max(logP); x0 = xTrain(:,x0);
  xEval = maxEnt_gibbs(nSamplesEval, burnIn, lambdaMPF, x0, modelFit);
 fxEval = setup_features_maxent(xEval', modelFit)';
 

 
 
%% depict results
%--------------------------------------------------------------------------
 disp('Visualizing results')
 figure(1); 
 subplot(1,2,1)
  if strcmp(modelTrue(1:5),'ising')
   plot(lambdaTrue, 'ko-', 'linewidth', 2); 
  end
  hold on  
  plot(lambdaMPF, 'bo-', 'linewidth', 2);
  line([d,d]+            0.5, [min([lambdaTrue;lambdaMPF]), max([lambdaTrue;lambdaMPF])], 'color', 'k')
  line([d,d]+(d*(d-1)/2)+0.5, [min([lambdaTrue;lambdaMPF]), max([lambdaTrue;lambdaMPF])], 'color', 'k')  
  hold off
  title('True (if existent) and recovered model parameters \lambda_i')
  box off; set(gca, 'TickDir', 'out')
  axis([0.5, size(fxEval,1), min([lambdaTrue;lambdaMPF]),max([lambdaTrue;lambdaMPF])])
  xlabel('i'); ylabel('\lambda_i')
  
 subplot(1,2,2)
  mfxTrain = mean(fxTrain(:, end/2:end),2); 
  plot(mfxTrain, 'ko-', 'linewidth', 2)
  hold on
  mfxEval = mean(fxEval(:,end/2:end),2);
  plot( mfxEval, 'go-', 'linewidth', 2)
  mfxEval = mean(fxEval,2);
  plot( mfxEval, 'bo-', 'linewidth', 2)
  line([d,d]+            0.5, [min([lambdaTrue;lambdaMPF]), max([lambdaTrue;lambdaMPF])], 'color', 'k')
  line([d,d]+(d*(d-1)/2)+0.5, [min([lambdaTrue;lambdaMPF]), max([lambdaTrue;lambdaMPF])], 'color', 'k')
  hold off
  title('True and recovered distribution moments E[f_i(X)]')
  box off; set(gca, 'TickDir', 'out')
  axis([0.5, size(fxEval,1), min([mfxTrain;mfxEval]),max([mfxTrain;mfxEval])])
  xlabel('i'); ylabel('E[f_i(X)]')  
  
  
 % check Markov chain 
 tmp = zeros(length(nSamplesTrain-1000)); % check for trends in count K
 for t = 1:nSamplesTrain-1000
  tmp(t) = sum(mean(xTrain(:,t:t+1000),2),1);
 end
 figure; plot(tmp); 