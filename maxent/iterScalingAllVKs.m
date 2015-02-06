function [delta_k, LL, deltaRecurse] = iterScalingAllVKs(Efx_K, Efy_K, VKold, delta0, fitoptions)

if nargin < 4 || isempty(delta0)
    delta0 = zeros(size(Efx_K));
end

if nargin < 5 || isempty(fitoptions)
    fitoptions = struct;
    fitoptions.optTol=1e-100; 
    fitoptions.progTol=1e-100; 
    fitoptions.display='off';
    fitoptions.MaxIter=3000;
    fitoptions.maxFunEvals=10000;
end

if ~isfield(fitoptions, 'maxK') || isempty(fitoptions.maxK)
 fitoptions.maxK = find(Efx_K==0, 1, 'first')-2; 
  % -1 for 'last consecutive Efx_K>0', -1 for K starting at 1, Efx_K at K=0  
  disp(['maxK = ',num2str(fitoptions.maxK)])
end
if ~isfield(fitoptions, 'sig2_l2') || isempty(fitoptions.sig2_l2)
 fitoptions.sig2_l2 = 10000;    
end
if ~isfield(fitoptions, 'sig2_sm') || isempty(fitoptions.sig2_sm)
 fitoptions.sig2_sm = 4 * max(var(VKold(2:fitoptions.maxK)), 0.001);   
end
if ~isfield(fitoptions, 'tau') || isempty(fitoptions.tau)
 fitoptions.tau = 3;   
end

idxBad = (Efy_K == 0); % finite sampling issue: if E_emp[f_k(X)] > 0 for 
 % the empirical distribution p_emp(X), but E[f_k(Y)] for the MCMC sample 
 % for any k, 0 <= k <= n, the most recent guess of lambda_true, 
 % then we are in trouble.
 % The equations for this would situation try to set V(K=k) to Infinity. 
 % We'll catch these cases. 
idxBad(fitoptions.maxK+2:end) = true; % second type of error catching: 
% Compute covariance matrix of Gaussian prior
Sigma = diag(fitoptions.sig2_l2 * ones(fitoptions.maxK,1)) + ...
        fitoptions.sig2_sm * ...
        exp(- 0.5*(ones(fitoptions.maxK,1)*(1:fitoptions.maxK) - ...
                  (1:fitoptions.maxK)'*ones(1,fitoptions.maxK)).^2/fitoptions.tau^2 );
Sigma0Rest = fitoptions.sig2_sm * exp(- 0.5*( (1:fitoptions.maxK).^2/fitoptions.tau^2 ) ); 
Sigma = Sigma - (Sigma0Rest'*Sigma0Rest)/(fitoptions.sig2_l2+fitoptions.sig2_sm);
SigmaVKold = Sigma \ VKold(1+(1:fitoptions.maxK));
% the above matrix is n-by-n (instead of n+1-by-n+1 for K = 0,...,N) 
% because delta_0 = delta_k(1) = 0 is fixed

[delta_kMax,LL,~,~] = minFunc( @LLVKs, delta0((1:fitoptions.maxK)+1), ...
                               fitoptions, ...
                               Efx_K(1:fitoptions.maxK+1), ...
                               Efy_K(1:fitoptions.maxK+1), ...
                               idxBad((1:fitoptions.maxK)+1), ...
                               Sigma, SigmaVKold);

delta_k = delta0; 
delta_k(1+(1:fitoptions.maxK)) = delta_kMax;

if nargout > 2
 deltaRecurse = sanityCheck(delta_k, Efx_K, Efy_K);
end

    
function [LL, dLL] = LLVKs(delta_k, Efx_K, Efy_K, idxBad, Sigma, SigmaVKold)
 % delta_k now contains delta-Terms for K = 1, ..., maxK. delta(K) for K=0
 % is again fixed to zero. delta(K) for K > maxK shall not be updated. 
 % idxBad tells which of the V(K)-Terms, K = 1,...,maxK are also not to be
 % updated. 
 % Efx_K, Efy_K still contain P(K=k)-Terms for K = 0, N !
  maxK = length(delta_k);
  Eexpdelta_k = Efy_K(1) + Efy_K(2:maxK+1)' * exp(delta_k);
  SigmaDelta = Sigma \ delta_k; % maxK-by-maxK * maxK-by1 = maxK-by-1 vector
  LL = - delta_k' * Efx_K(2:maxK+1) + ...
         log(Eexpdelta_k) + ...
         delta_k' * (0.5*SigmaDelta + SigmaVKold);
 dLL = - Efx_K(2:maxK+1) + ...
        (Efy_K(2:maxK+1) .* exp(delta_k)) / Eexpdelta_k + ...
         SigmaDelta + SigmaVKold;
 dLL(idxBad) = 0; % do not update these
end

function res = sanityCheck(delta_k, Efx_K, Efy_K)
  res = log(Efx_K ./ Efy_K) + log(Efy_K' * exp(delta_k));
end

end