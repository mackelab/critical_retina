function MSE = getMSEtraceOldRun(fname, n)


load('/home/marcel/criticalityIterScaling/data/EfxSimData.mat')
Efx = Efx{n}; % get empirical E[f(X)] for all subpopulations of size 10*n

n = 10*n; % actual population size 

fname = [fname, 'n', num2str(n)];
cd(['/home/marcel/criticalityIterScaling/'])

if exist([fname, 'asfinal.mat'], 'file')
  disp('loading file')
  load([fname, 'asfinal.mat']);
   % gives lambdaHat, fD (fit diagnostics), n, EfxHat, EE, varE
else % assemble data from currently most recent iteration step
  warning('run results not found')  
end


fDescrJ = nchoosek(1:n,2)'; 

idxRun = zeros(length(fD),1);
maxIter = 0;
for i = 1:length(idxRun)
if ~isempty(fD{i})
 idxRun(i) = 1;
 maxIter = max([maxIter, size(fD{i}.EfyTrace,2)]);% maximum number of iters
end                                              % for ALL data sets

end                                               
idxRun = find(idxRun);
maxRun = max(idxRun); 

MSE.h=Inf*ones(maxIter,maxRun); MSE.V=Inf*ones(maxIter,maxRun);MSE.cov=Inf*ones(maxIter,maxRun);
for idx = 1:length(idxRun) % go through data sets
 j = idxRun(idx); % data set index (in case it differs from idx)
 
 maxIter=size(fD{j}.EfyTrace,2); %maximum number of iters for THIS data set

 Efy = fD{j}.EfyTrace; 
 covs.x = Efx((n+1):(n*(n+1)/2),j) - ...
         (Efx(fDescrJ(1, :),j).* Efx(fDescrJ(2, :),j));
 for i = 1:maxIter
% Compute errors on E[f(X)], split into blocks corresponding to (h,J,V)
    MSE.h(i,j) = mean( (Efx(1:n,j)-Efy(1:n,i)).^2 );
    MSE.V(i,j) = mean( (Efx(end-n:end,j)-Efy(end-n:end,i)).^2 ); 
    MSE.J(i,j) = mean( (Efx((n+1):(n*(n+1)/2),j)-Efy((n+1):(n*(n+1)/2),i)).^2 ); 
    
    covs.y = Efy((n+1):(n*(n+1)/2),i) - ...              % cov(a,b) = E(a*b) -
          (Efy(fDescrJ(1, :),i).* Efy(fDescrJ(2, :),i)); %          E(a)*E(b)
    MSE.cov(i,j) = mean( (covs.x - covs.y).^2 );
 end
end

end % end function
