function [lambdaHat, fD] = iterScaling(xTrain, fitoptions, beta, eps, fname, ifSave, hJV, ifbwVK)
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
%      .sig2_l2: regularization strength for the ridge-part of the bwVK
%      .sig2_sm: regularization strength for the smoothing-part of the bwVK
%   fname: (optional) string of file name for intermediate result storage
%   ifSave: boolean specifying whether or not to save intermediate results
%   hJV: 3-by-1 vector of booleans giving whether to include h, J and/or V 
%   ifbwVK: boolean specifying whether or not to update V(K) block-wise
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

ticTime = now;

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

if nargin < 7 || isempty(hJV)
  hJV = ones(3,1); % default: include terms for h, J and V(K) (full model)
end
if nargin < 6 || isempty(ifSave)
  ifSave = true; % default: save intermediate results after each iteration
end
if nargin < 5 % note that setting fname = [] is possible for test purposes
  fname = date;  % default: save with current date
  ifSave = false;
end
if nargin < 4 || isempty(eps)
  eps = 10^(-2) * ones(3,1); % default error tolerance: 1%
end
if nargin < 3 || isempty(beta)
  beta = zeros(size(Efx)); % default: no regularization
end

if all(size(fitoptions.nSamples)<2) % if nSamples is constant for all iters
  fitoptions.nSamples = fitoptions.nSamples*[0;ones(fitoptions.maxIter,1)];
end
if all(size(fitoptions.burnIn)<2)   % if burnIn is constant for all iters
  fitoptions.burnIn = fitoptions.burnIn * [0;ones(fitoptions.maxIter,1)];
end

% fitoptions for the call to minFunc during iterScalingAllVKs()
fitoptionsbwVK.maxK = find(Efx(end-n:end)==0, 1, 'first')-2;
fitoptionsbwVK.sig2_l2 = 400;   
fitoptionsbwVK.sig2_sm = []; % will set this for each individual call   
fitoptionsbwVK.tau = 2;   
fitoptionsbwVK.optTol=1e-100; 
fitoptionsbwVK.progTol=1e-100; 
fitoptionsbwVK.display='off';
fitoptionsbwVK.MaxIter=3000;
fitoptionsbwVK.maxFunEvals=10000;

idxBad = (Efx == 0); % identify pathological cases where lambda_i -> -Inf

if ~hJV(1) % declare parameters of h to be 'bad', i.e. do not fit them
 idxBad = idxBad | ((1:length(Efx)) <= n)';  
 fitoptions.lambda0(1:n,:) = 0;
 eps(1) = Inf; % ignore errors for h-Terms when checking for convergence
end
if ~hJV(2) % declare parameters of J to be 'bad', i.e. do not fit them
 idxBad = idxBad | (((1:length(Efx))>n)&((1:length(Efx))<length(Efx)-n))';  
 fitoptions.lambda0((n+1):(n*(n+1)/2),:) = 0; 
 eps(2) = Inf; % ignore errors for J-Terms when checking for convergence
end
if ~hJV(3) % declare parameters of V to be 'bad', i.e. do not fit them
 idxBad = idxBad | ((1:length(Efx)) >= length(Efx)-n)';  
 fitoptions.lambda0((n*(n+1)/2+1):end,:) = 0; 
 eps(3) = Inf; % ignore errors for V-Terms when checking for convergence
end

%fDescr  = [[1:n; (1:n)*nan],nchoosek(1:n,2)',[0:n; (0:n)*nan]]; 
fDescrJ = nchoosek(1:n,2)'; 

covs.x = Efx((n+1):(n*(n+1)/2)) - ...              % covariance as computed
        (Efx(fDescrJ(1, :)).* Efx(fDescrJ(2, :))); % from Efx

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
    disp('No inizialitaion for \lambda given. Generating new initialization')
    lambdaHat(:,1) = [randn(n,1);randn(n*(n-1)/2,1)/sqrt(n);randn(n+1,1)];            
  else
    disp('Using user-provided initialization for \lambda')
    lambdaHat(:,1) = fitoptions.lambda0;                                                 
  end
  
  lambdaHat(idxBad,1) = 0;     % Set those parameters for features f_i with
  lambdaHat(n*(n+1)/2+1,1) = 0;% E[f_i]=0 and the feature for K=0 to zero 

  lambdaHat(end-n+fitoptionsbwVK.maxK+1:end,1) = -1000; 
  
  % Generate first MCMC chain element x0 using E[X] from a maxEnt model 
  EX = exp(lambdaHat(1:n,1))./(1+exp(lambdaHat(1:n,1))); % with only h,
  x0(:,1) = double(rand(n,1)<EX);              % i.e. no parameters J, L
  
  fD.idxBad = idxBad; % indexes of parameters not to be touched
  fD.Efx = Efx; % what we try to achieve
  
  minIter = 2; % standard starting index of MCMC chain (iter=1 is x0)

if ~isempty(fname)
  % Check for previous results from potentially cancelled run
  cd(['/home/marcel/criticalityIterScaling/results/',fname,'/'])
  names = strsplit(ls); % may take forever...
  names = sort(names);
  if ~isempty(names) && ~strcmp(names{end},'') % '' is always first after sorting
   lfile = names{end-1}; % load second-last iteration (last file could be corrupted)
   minIter = str2num(lfile(end-8:end-4)); 
   if minIter < fitoptions.maxIter
    load(lfile); % 
    cd '/home/marcel/criticalityIterScaling/' % get out of crowded /results!
     % gives Efy, deltaIter, deltaLL, idxIter, lambdaIter, x0Iter
    x0(:,minIter-1) = x0Iter; 
    lambdaHat(:,minIter-1) = lambdaIter; 
    clear Efy deltaIter deltaLL idxIter lambdaIter x0Iter
    disp(['restarting from earlier run ', fname, ', iter #', num2str(minIter), '/', num2str(fitoptions.maxIter)])
   end
  end
end   

% MAIN LOOP

  for iter = minIter:fitoptions.maxIter+1
      
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
          
          if ifbwVK
             fitoptionsbwVK.sig2_sm = ... % get strength of V(K)-smoothing 
                   var(lambdaHat(end-n+1:end-n+fitoptionsbwVK.maxK,iter));               
             [deltaVK, deltaLLVK] = iterScalingAllVKs(Efx(end-n:end),...
                                                      Efy(end-n:end),...
                                           lambdaHat(end-n:end,iter),...
                                                      zeros(n+1,1),  ...
                                                      fitoptionsbwVK);
              delta(end-n:end) = deltaVK;
          end
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

          if ifbwVK
             [deltaVK, deltaLLVK] = iterScalingAllVKs(Efx(end-n:end),...
                                                      Efy(end-n:end),...
                                           lambdaHat(end-n:end,iter),...
                                                      zeros(n+1,1),  ...
                                                      fitoptionsbwVK);

              delta(end-n:end) = deltaVK;
          end
          
        otherwise   % do nothing as for no regularization, but warn       
          delta = log( (Efx .* ( 1 - Efy )) ./ (Efy .* ( 1 - Efx )) );            
          if ifbwVK
             [deltaVK, deltaLLVK] = iterScalingAllVKs(Efx(end-n:end),...
                                                      Efy(end-n:end));
              delta(end-n:end) = deltaVK;
          end          
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
        case 'l1'   % add l1 regularization terms where applicable
          deltaLL(1:(n*(n+1)/2)) = deltaLL(1:(n*(n+1)/2)) + ... 
           beta(1:(n*(n+1)/2)) .* ...
               (abs(lambdaHat(1:(n*(n+1)/2),iter)+delta(1:(n*(n+1)/2))) ...
                 - abs(lambdaHat(1:(n*(n+1)/2),iter))); 
          if ifbwVK  % l2-regularized block-wise update of V(K)
            deltaLL(end-n:end) = deltaLLVK;             
          else % if not l2-regularizing V(K), then treat as all others
            deltaLL(end-n:end) = deltaLL(end-n:end) + ... % l1-regularize
             beta(end-n:end) .* (abs(lambdaHat(end-n:end,iter)+delta) ...
                                 - abs(lambdaHat(end-n:end,iter))); 
          end            
        otherwise   % do nothing as for no regularization, but warn          
          disp('Unknown regularization chosen. Solution not regularized')
    end
    deltaLL(idxBad)      = Inf; % never update 'bad' components of lambda
    deltaLL(n*(n+1)/2+1) = Inf; % or the weight of the feature for K=0
    
    % Pick candidate update that gives highest gain
    [~, idxj(iter)] = min(deltaLL); % minimizing NEGATIVE log-likelihood  
    deltaLLs(iter-1) = deltaLL(idxj(iter)); % store for debugging purposes
    
    
    % Update correct component of parameter vector lambda
    if ifbwVK && idxj(iter) > n*(n+1)/2 % update V(K) as single block
      lambdaHat(end-n:end,iter)  = lambdaHat(end-n:end,iter) + deltaVK;
    else % either no blocked V(K)-update or updating a h- or J-term
      lambdaHat(idxj(iter),iter) = lambdaHat(idxj(iter),iter) + ...
                                                         delta(idxj(iter));         
    end
       
    % Compute errors on E[f(X)], split into blocks corresponding to (h,J,V)
    MSE.h = mean( (Efx(1:n)-Efy(1:n)).^2 );
    MSE.V = mean( (Efx(end-n:end)-Efy(end-n:end)).^2 ); 
    MSE.J = mean( (Efx((n+1):(n*(n+1)/2))-Efy((n+1):(n*(n+1)/2))).^2 ); 
    
    covs.y = Efy((n+1):(n*(n+1)/2)) - ...              % cov(a,b) = E(a*b) -
          (Efy(fDescrJ(1, :)).* Efy(fDescrJ(2, :))); %            E(a)*E(b)
    MSE.cov = mean( (covs.x - covs.y).^2 );

    MSEperc.h = sqrt(MSE.h)/sqrt(mean(Efx(1:n).^2));        % percentage of
    MSEperc.V = sqrt(MSE.V)/sqrt(mean(Efx(end-n:end).^2));  % mean-squares
    MSEperc.cov = sqrt(MSE.cov)/sqrt(mean(abs(covs.x).^2)); % of Efx
    
    % Save
    if ifSave
     
     idxIter = idxj(iter);
     x0Iter  = x0(:,iter);
     lambdaIter = lambdaHat(:,iter);
     deltaIter = delta;
     fnames = [fname,'_Iter_',num2str(iter, '%0.5d')];
     fnames=['/home/marcel/criticalityIterScaling/results/',fname,...
             '/', fnames,'.mat'];
     save(fnames, 'deltaLL', 'deltaIter', 'idxIter', 'Efy', 'x0Iter', ...
                             'lambdaIter', 'MSE', 'MSEperc', 'covs') 
    end
    
    % Check for convergence
    % 1. Ran for more than 24 h hours
    tocTime = now; 
    if (tocTime - ticTime) > 1
      break;
    end
    % 2. Did more updates than 10 times the number of parameters
    if iter > (10 * n * (n+3) / 2)
      break;
    end
    % 3. Error smaller than tolerance    
    if ( (MSEperc.h < eps(1)) && ...   % check h, V first, 
         (MSEperc.V < eps(3)) && ... % as is faster
         (MSEperc.cov  < eps(2)) ) 
          break;
    end
   
    
  end % END MAIN LOOP

     fD.deltaLL = deltaLL;    % currently possible gains in log-likelihood
     fD.deltaLLs = deltaLLs(1:iter-1);  % trace of realized gains in log-likelihood
     fD.lambdaTrace = lambdaHat(:,2:iter); % trace of parameter estimates
     fD.idxUpdate = idxj(2:iter); % trace of parameters picked for updating
     fD.deltas = deltas(:,1:iter-1);  % trace of sizes of changes in parameters
     fD.EfyTrace = Efys(:,2:iter);  % trace of resulting expected values
     fD.Efy = Efy; % what we did achieve in quality up to this iteration
     fD.x0 = x0(:,1:iter); % trace of initial chain elements for each MCMC draw
     fD.MSE = MSE;         % mean-squared errors on final results, 
     fD.MSEperc = MSEperc; % absolute and relative to E[f(X)] magnitudes
     fD.maxK = fitoptionsbwVK.maxK;
     fD.eps = eps;
     fD.fitoptions = fitoptions;
     fD.fitoptionsbwVK = fitoptionsbwVK;
     fD.nIters = iter-1; % 'first iteration' is initialization
end

lambdaHat = lambdaHat(:,fitoptions.maxIter+1);

end