function [lambdaHat, fitDiagnostics] = iterScaling(xTrain, fitoptions)
% code for fitting a maxEnt model with (improved?) iterative scaling

% 1. Input formatting and pre-computations
%------------------------------------------
[n, ~] = size(xTrain);    
if size(xTrain,2) > 1 % if feeding training data directly
 fxTrain = setup_features_maxent(xTrain', fitoptions.modelFit);
 Efx = full(mean(fxTrain,1)');
else                  % if feeding expectation of features E[f(X)] directly
 n = round(sqrt(2*n + 0.25) - 1.5); % #dimensions to #features, backwards
 Efx = xTrain;
end
clear fxTrain % may take a whole lot of memory 

idxBad = (Efx == 0);

if ~isempty(fitoptions.lambda0) && size(fitoptions.lambda0,2) ~= fitoptions.nRestart
  fitoptions.lambda0=repmat(fitoptions.lambda0(:,1),1,fitoptions.nRestart);
end

x0 = n; 

%Efx(n*(n+1)/2+1) = [];                                                                   %TEST
%fitoptions.lambda0(n*(n+1)/2+1) = [];                                                    %TEST

%% 2. Start the update iteration
%------------------------------------------
lambdaHat = zeros(length(Efx),fitoptions.maxIter+1);
d = zeros(fitoptions.maxIter,1);
idxj = zeros(fitoptions.maxIter,1);
deltas = zeros(size(lambdaHat));
deltaLLs = zeros(fitoptions.maxIter,1);
for r = 1:fitoptions.nRestart
  if isempty(fitoptions.lambda0) || size(fitoptions.lambda0,1) ~= size(lambdaHat,1)
    lambdaHat(:,1) = [randn(n,1);randn(n*(n-1)/2,1)/sqrt(n);randn(n+1,1)];                %TEST
  else
    lambdaHat(:,1) = fitoptions.lambda0;                                                 
  end
  
  lambdaHat(idxBad,1) = 0;    % Set those parameters for features f_i with
  lambdaHat(n*(n+1)/2+1) = 0; % E[f_i]=0 and the feature for K=0 to zero 
    
  for iter = 2:fitoptions.maxIter+1
   disp([num2str(iter-1), '/' num2str(fitoptions.maxIter)])
    lambdaHat(:,iter)  = lambdaHat(:,iter-1);
    [Efy,~,x0] = maxEnt_gibbs_pair_C(fitoptions.nSamples, fitoptions.burnIn, lambdaHat(:,iter), x0, fitoptions.machine);
    
    %for innerIter = 1:fitoptions.maxInnerIter
    delta = log( (Efx .* ( 1 - Efy )) ./ (Efy .* ( 1 - Efx )) );
    deltas(:,iter) = delta;

    
    deltaLL= - delta .* Efx + log( 1 + (exp(delta)-1) .* Efy ) ; 
    deltaLL(idxBad)      = Inf; % do not update 'bad' components of lambda
    deltaLL(n*(n+1)/2+1) = Inf;
    
    [~, idxj(iter)] = min(deltaLL);   
    
    lambdaHat(idxj(iter),iter) = lambdaHat(idxj(iter),iter) + delta(idxj(iter)); 
    
    %end
  end
%  [Efy,~,x0] = maxEnt_gibbs_pair_C(10*fitoptions.nSamples, fitoptions.burnIn, lambdaHat(:,fitoptions.maxIter), x0, fitoptions.machine);
%   delta = log( (Efx .* ( 1 - Efy )) ./ (Efy .* ( 1 - Efx )) );
%   idxGood = delta<0 & ~idxBad;
%   idxGood(n*(n+1)/2+1) = false;
%   delta(~idxGood) = 0;
%   deltas(:,end) = delta;
%   lambdaHat(:,end) = lambdaHat(:,fitoptions.maxIter) + delta; 
end

fitDiagnostics.deltaLL = deltaLL;
fitDiagnostics.lambdaTrace = lambdaHat(:,2:end); 
fitDiagnostics.idxUpdate = idxj(2:end);
fitDiagnostics.sizeUpdate = d; 
fitDiagnostics.deltas = deltas;
fitDiagnostics.idxBad = idxBad;
fitDiagnostics.Efx = Efx; % what we tried to achieve
fitDiagnostics.Efy = Efy; % what we did achieve
lambdaHat = lambdaHat(:,fitoptions.maxIter);

% figure; 
% subplot(2,3,1),
% plot(lambdaHat(:,:)); hold on; plot(lambdaHat(:,end), 'r', 'linewidth', 3); plot(lambdaTrue,'k', 'linewidth', 2)
% subplot(2,3,2),
% plot(lambdaHat(:,:)'); hold on; plot(1.1*fitoptions.maxIter*ones(size(lambdaTrue)), lambdaTrue,'k*', 'markerSize', 5)
% subplot(2,3,3);
% plot(Efx, 'b*-'); hold on; plot(Efy, 'co-');
% subplot(2,3,4),
% plot(d); box off; 
% subplot(2,3,5)
% h = histc(idxj, 0:length(Efx)); bar(0:length(Efx), h/sum(h));

end