function assembleIterScalingQuick(fname, n, idxRepet)

d = n; %variable conventions...

fD = cell(length(idxRepet),1);
lambdaHat = zeros(n*(n+3)/2+1, length(idxRepet));

for i = idxRepet
    
  cd '/home/marcel/criticalityIterScaling/results/'
  names = strsplit(ls([fname, 'idxRepet', num2str(i),'*'])); % may take forever...
  cd '/home/marcel/criticalityIterScaling/'
 
  names = sort(names);
  iteri = zeros(length(names),1);
  for j = 2:lenth(names)
    iteri(j) = str2num(names{j}(end-8:end-4));
  end
  
  nf = d*(d+3)/2+1; % number of featurs for d = 100
  cmi = max(iteri); % current maximum iteration count
     
  fD{i}.deltaLLs = zeros(nf,cmi);  % trace of possible gains in log-likelihood
  fD{i}.deltas = zeros(nf,cmi);    % trace of sizes of changes in parameters
  fD{i}.lambdaTrace=zeros(nf,cmi); % trace of parameter estimates
  fD{i}.EfyTrace = zeros(nf,cmi);  % trace of resulting expected values
  fD{i}.x0 = zeros(d, cmi);        % trace of initial chain elements 
  fD{i}.idxUpdate = zeros(1,cmi);  % trace of parameters picked for updating
  fD{i}.deltaLL = zeros(1,cmi);    % trace of realized gains in log-likelihood
  
   disp('loading data and attaching')
  for j = 2:length(iteri)
     idx = iteri(j); % idx may differ from j if some files are missing
      disp([ 'loading data ', num2str(j), '/', num2str(length(iteri))] )
     load(['/home/marcel/criticalityIterScaling/results/',names{i}]);         % now load stuff and collect!
     fD{i}.deltaLLs(:,idx)    = deltaLL;
     fD{i}.deltas(:,idx)      = deltaIter;
     fD{i}.lambdaTrace(:,idx) = lambdaIter;
     fD{i}.EfyTrace(:,idx)    = Efy;
     fD{i}.x0(:,idx)          = x0Iter;
     fD{i}.idxUpdate(idx)     = idxIter;
     fD{i}.deltaLL(idx)       = deltaLL(idxIter);
  end
  fD{i}.Efy = Efy;   % what we did achieve in quality by the end
  lambdaHat{i} = fD.lambdaTrace(:,idx);

end

  disp('Storing results up to n-ith step')     
     save([fname, 'final.mat'], 'lambdaHat', 'fD')         % actual outcome
       