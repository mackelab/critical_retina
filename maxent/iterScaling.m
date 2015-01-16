function [lambdaHat, fD] = iterScaling(xTrain, fitoptions, beta, fname)
% code for fitting a maxEnt model with (improved?) iterative scaling.
% Minimizes the negative log-likelihood of the data in xTrain under the 
% maximum entropy model by (blockwise) coordinate descent. 
% 
% input:
%   xTrain:     n-by-N matrix of N many n-dimensional training patterns x 
%               OR #features-by-1 vector of sample means E_emp[f(X)]
%   fitoptions: structure containing options for the fitting procedure
%      .modelFit:  string that specifies which model was used
%      .regular:   string that specifies which regularization is to beused
%      .lambda0:   #features-by-1 vector of initial values for lambda
%      .nRestarts: number of complete algorithm restarts (test convexity)
%      .maxIter:   maximal number of allowed main loop iterations
%      .maxInnerIter: maximum number of allowed inner loop iterations
%      .nSamples:  number of Gibbs samples for each new MCMC sample
%      .burnIn  :  number of Gibbs samples to be discarded from MCMC sample
%      .machine :  string specifying which Gibbs implementation is used
%   fname: (optional) string of file name for intermediate result storage
% output:
%   lamdaHat:   #features-by-1 vector of lambda as estimated from data
%   fD:         structure containing several fitting diagnostics 
%      .deltaLL:     full sequence of changes in negative log-likelihood
%      .idxBad:      list of 'bad' feature dimensions (where E[f_i(X)] = 0)
%      .lambdaTrace: all intermediate estimates of lambda
%      .idxUpdate:   sequence of dimensions of lambda updated during fit
%      .deltas:      all candidate optimal update step sizes
%      .Efx:         sought-after data means E_emp[f(X)]
%      .Efy:         actually returned model means E_lambda[f(X)] 


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

if nargin < 4
  fname = date;
end
if nargin < 3
  beta = zeros(size(Efx));
end

if any(beta >= Efx | beta >= 1-Efx)
  %beta(beta>Efx) = 0.9* min([1-Efx(beta>Efx),Efx(beta>Efx)],[],2); 
  disp('Warning: Some of the regularizers \beta_j were chosen larger than E[f_j(X)]!')
end

if all(size(fitoptions.nSamples)<2)
  fitoptions.nSamples = fitoptions.nSamples*[0;ones(fitoptions.maxIter,1)];
end
if all(size(fitoptions.burnIn)<2)
  fitoptions.burnIn = fitoptions.burnIn * [0;ones(fitoptions.maxIter,1)];
end
%Efxh = Efx(1:n);               % appears wise to split the feature
%EfxJ = Efx((n+1):(n*(n+1)/2)); % dimensions of f(X) up according to the
%EfxV = Efx((end-n):end);       % upcoming block of coordindate descent

if ~isempty(fitoptions.lambda0) && size(fitoptions.lambda0,2) ~= fitoptions.nRestart
  fitoptions.lambda0=repmat(fitoptions.lambda0(:,1),1,fitoptions.nRestart);
end % copy initialization for every run if not already provided

idxBad = (Efx == 0); % identify pathological cases where lambda_i -> -Inf



                
%% 2. Start the update iteration
%------------------------------------------
lambdaHat = zeros(length(Efx),fitoptions.maxIter+1);
Efys      = zeros(size(lambdaHat)); % keeps track of expectation goodness
x0   = zeros(n, fitoptions.maxIter+1); % keeps track of MCMC chain starts
idxj = zeros(fitoptions.maxIter,1); % keeps track of updated dimensions
deltas = zeros(size(lambdaHat));    % keeps track of optimal step sizes
deltaLLs = zeros(fitoptions.maxIter,1); % keeps track of LL improvements




for r = 1:fitoptions.nRestart
    
  % generate / load parameter initialization
  if isempty(fitoptions.lambda0) || size(fitoptions.lambda0,1) ~= size(lambdaHat,1)
    lambdaHat(:,1) = [randn(n,1);randn(n*(n-1)/2,1)/sqrt(n);randn(n+1,1)];            
  else
    lambdaHat(:,1) = fitoptions.lambda0;                                                 
  end
  
  lambdaHat(idxBad,1) = 0;     % Set those parameters for features f_i with
  lambdaHat(n*(n+1)/2+1,1) = 0;% E[f_i]=0 and the feature for K=0 to zero 

  % Generate first MCMC chain element x0 using E[X] from a maxEnt model 
  EX = exp(lambdaHat(1:n,1))./(1+exp(lambdaHat(1:n,1))); % with only h,
  x0(:,1) = double(rand(n,1)<EX);              % i.e. no parameters J, L
  
  fD.idxBad = idxBad; % indexes of parameters not to be touched
  fD.lambdaInit = lambdaHat(:,1);
  fD.Efx = Efx; % what we try to achieve
  
% MAIN LOOP

  for iter = 2:fitoptions.maxIter+1
      
    disp([num2str(iter-1), '/' num2str(fitoptions.maxIter)])
   
    lambdaHat(:,iter)  = lambdaHat(:,iter-1);
    [Efy,~,x0(:,iter)] = maxEnt_gibbs_pair_C(fitoptions.nSamples(iter), ...
                                             fitoptions.burnIn(iter), ...
                          lambdaHat(:,iter), x0(:,iter-1), fitoptions.machine);
    % returning x0 (now a full n-by-1 vector) for use as initial element of
    % the next Gibbs chain, assuming that the distributions do not change 
    % much from one iteration to the next. 
    Efys(:,iter) = Efy;
    %Efyh = Efy(1:n);               % appears wise to split the feature
    %EfyJ = Efy((n+1):(n*(n+1)/2)); % dimensions of f(Y) up according to the
    %EfyV = Efy((end-n):end);       % upcoming block of coordindate descent
    
    %for innerIter = 1:fitoptions.maxInnerIter % using sample Y several times
    
    % Compute optimal candidate step lengths for each dimension of lambda
    switch fitoptions.regular
        case 'none' % immediately done
          delta = log( (Efx .* ( 1 - Efy )) ./ (Efy .* ( 1 - Efx )) );
        case 'l1'
          delta0 = log( (1-Efy) ./ Efy ); % shows up twice below,
                                          % so compute once and store
          % \delta_+ part                                
          delta =      delta0     + log( (Efx-beta)./(1 - Efx + beta) );
          % find out where using \delta_+ is NOT appropriate
          ic = (beta >= Efx) | (delta <= -lambdaHat(:,iter)); 
          % \delta_- part
          delta(ic) =  delta0(ic) + log( (Efx(ic)+beta(ic))./(1-Efx(ic)-beta(ic)));
          % find out where using \delta_+ OR \delta_- is NOT appropriate
          ic = (ic & ((beta >= 1-Efx) | (delta >= -lambdaHat(:,iter))) );
          % \delta_0 part
          delta(ic) = -lambdaHat(ic,iter);
          
        otherwise   % do nothing as for no regularization, but warn       
          delta = log( (Efx .* ( 1 - Efy )) ./ (Efy .* ( 1 - Efx )) );            
          disp('Unknown regularization chosen. Solution not regularized')
    end

    
    deltas(:,iter-1) = delta;
    
    % Preliminary code for block-wise updates
    %deltah = log( (Efxh ./ ( 1 - Efxh )) .* ...
    %              (1-Efyh+(exp(deltah(j)))) );
    %deltaJ = log( (EfxJ .* ( 1 - EfyJ )) ./ (EfyJ .* ( 1 - EfxJ )) );
    %deltaV = 0;                                                                      % FILL IN HERE
    %deltas(1:n) = deltah;                      % store 
    %deltas((n+1):(n*(n+1)/2),iter-1) = deltaJ; % for debugging 
    %deltas((end-n):end) = deltaV;              % purposes
    
    % Compute gains in log-likelihood for these candidate updates
    deltaLL = - delta .* Efx + log( 1 + (exp(delta)-1) .* Efy ) ;     
    switch fitoptions.regular
        case 'none' % immediately done
        case 'l1'
          deltaLL = deltaLL + ... % add regularization terms
           beta .* (abs(lambdaHat(:,iter)+delta) - abs(lambdaHat(:,iter))); 
        otherwise   % do nothing as for no regularization, but warn          
          disp('Unknown regularization chosen. Solution not regularized')
    end
    deltaLL(idxBad)      = Inf; % never update 'bad' components of lambda
    deltaLL(n*(n+1)/2+1) = Inf; % or the weight of the feature for K=0
    
    % Pick candidate update that gives highest gain
    [~, idxj(iter)] = min(deltaLL); % minimizing NEGATIVE log-likelihood  
    deltaLLs(iter-1) = deltaLL(idxj(iter)); % store for debugging purposes
    
    
    % Update correct component of parameter vector lambda
    lambdaHat(idxj(iter),iter) = lambdaHat(idxj(iter),iter) + delta(idxj(iter)); 

    % Save
     fD.deltaLL = deltaLL;    % currently possible gains in log-likelihood
     fD.deltaLLs = deltaLLs;  % trace of realized gains in log-likelihood
     fD.lambdaTrace = lambdaHat(:,2:end); % trace of parameter estimates
     fD.idxUpdate = idxj(2:end); % trace of parameters picked for updating
     fD.deltas = deltas;      % trace of sizes of changes in parameters
     fD.EfyTrace = Efys(:,2:end);  % trace of resulting expected values
     fD.Efy = Efy; % what we did achieve in quality up to this iteration
     fD.x0 = x0;   % trace of initial chain elements for each MCMC draw
     
     idxIter = idxj(iter);
     x0Iter  = x0(:,iter);
     lambdaIter = lambdaHat(:,iter);
     deltaIter = deltas(:,iter);
    fnames = [fname,'_Iter_',num2str(iter, '%0.5d')];
    fnames=['/home/marcel/criticalityIterScaling/results/',fnames,'.mat'];
    save(fnames, 'deltaLL', 'deltaIter', 'idxIter', 'Efy', 'x0Iter', 'lambdaIter') 
     
%     if ismember(idxj(iter), find(ic))
%         disp(['l1 set weight ', num2str(idxj(iter)), 'to 0!'])
%         disp(num2str(lambdaHat(idxj(iter),iter)))
%         disp(num2str(delta(idxj(iter))))
%         disp('--- deltaLL:')
%         disp(num2str(deltaLL))
%         disp('--- delta: ')
%         disp(num2str(delta))
%         disp('--- Efy: ')
%         disp(num2str(Efy))
%     end

    %end
    
  end % END MAIN LOOP
  
% Final step to enforce updating all components of lambda:  
%  [Efy,~,x0] = maxEnt_gibbs_pair_C(10*fitoptions.nSamples, fitoptions.burnIn, lambdaHat(:,fitoptions.maxIter), x0, fitoptions.machine);
%   delta = log( (Efx .* ( 1 - Efy )) ./ (Efy .* ( 1 - Efx )) );
%   idxGood = delta<0 & ~idxBad;
%   idxGood(n*(n+1)/2+1) = false;
%   delta(~idxGood) = 0;
%   deltas(:,end) = delta;
%   lambdaHat(:,end) = lambdaHat(:,fitoptions.maxIter) + delta; 
end

lambdaHat = lambdaHat(:,fitoptions.maxIter+1);

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