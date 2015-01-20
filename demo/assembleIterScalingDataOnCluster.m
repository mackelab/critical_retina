function assembleIterScalingDataOnCluster(runIDs)

%input: 
% - runIDs: subset of run IDs (first run was sessions numbered 1 to 8,
%           second run was session number 91 to 98) to assemble data from.

% IMPORTANT: don't try to run this on session groups that have different
%            numbers of digits assigned to them, e.g. keep apart 
%            sessions 1 to 8 from sessions 91 to 98 (stupid, but necessary
%            given the simple way of identifying session names below)
%            runIDs=[1,6:8] or runIDs=[91:98] is OK, runIDs=[1,91] is NOT!

disp('moving into data folder - please press key')
pause; 
cd('/home/marcel/criticalityIterScaling/results/')
%cd('C:\Users\Loki\Desktop\Criticality\Results\ClusterRuns\')
d = 100;
disp('making ls list')
dir = ls; % list of all files in results folder

disp('moving out of data folder')
cd('/home/marcel/criticalityIterScaling/') % getting out of that folder as
                                           % Matlab seems to be really mad
                                           % about all the files in there

disp('splitting ls list')
dirTmp = strsplit(dir);
m = 1;
disp('obtaining width of ls list')
for i = 1:length(dirTmp)
 m = max([m, length(dirTmp{i})]);
end
disp('reshaping ls list into long list (n-by-m) of filenames')
dir = reshape( blanks(m*length(dirTmp)), [length(dirTmp),m]);
for i = 1:length(dirTmp)
 dir(i,1:length(dirTmp{i})) = dirTmp{i};
end
    
range = 1:(4+length(num2str(runIDs(1)))); 
%dirStarts = reshape( blanks(length(range)*size(dir,2)), [size(dir,2),length(range)]);
dirStarts = cell(size(dir,1),1); % extract first letters of each file 
disp('obtaining initial parts of list names')
for i = 1:size(dir,1)            % to extract run ID 
     dirStarts{i} = dir( i, range );   
end

for i = runIDs % e.g. 1:8   % for each stored run ...
    
 disp(['obtaining numbers for session sess', num2str(i)])
  tStart = ['sess', num2str(i)];  % assumed first part of run idx
  diri = dir( strcmp(tStart, dirStarts), : ); % all save file names
  runIdx = diri(end,1:find(diri(end,:)=='.')-1); % full run idx (with nSamples)
  diri = diri(1:end-1,:); % first file contained pars, fitoptions etc.
  iteri = zeros(size(diri,1), 1); % storage for iteration counts

 disp('converting strings to numbers')
  for j = 1:size(diri,1)    % ... for each stored iteration within i'th run
    last_ = find(diri(j,:)=='.', 1, 'last'); % only one full stop per fname
    iteri(j) = str2double(diri(j,last_-5:last_-1)); % iteration count
  end
  nf = d*(d+3)/2+1; % number of featurs for d = 100
  cmi = max(iteri); % current maximum iteration count
     
  fD.deltaLLs = zeros(nf,cmi);  % trace of possible gains in log-likelihood
  fD.deltas = zeros(nf,cmi);    % trace of sizes of changes in parameters
  fD.lambdaTrace=zeros(nf,cmi); % trace of parameter estimates
  fD.EfyTrace = zeros(nf,cmi);  % trace of resulting expected values
  fD.x0 = zeros(d, cmi);        % trace of initial chain elements 
  fD.idxUpdate = zeros(1,cmi);  % trace of parameters picked for updating
  fD.deltaLL = zeros(1,cmi);    % trace of realized gains in log-likelihood
     
 disp('loading data and attaching')
  for j = 1:size(diri,1)-1
     idx = iteri(j); % idx may differ from j if some files are missing
     fname = diri(j,:); fname(fname==' ') = [];
      disp([ 'loading data ', num2str(j), '/', num2str(size(diri,1))] )
     load(['/home/marcel/criticalityIterScaling/results/',fname]);         % now load stuff and collect!
     fD.deltaLLs(:,idx)    = deltaLL;
     fD.deltas(:,idx)      = deltaIter;
     fD.lambdaTrace(:,idx) = lambdaIter;
     fD.EfyTrace(:,idx)    = Efy;
     fD.x0(:,idx)          = x0Iter;
     fD.idxUpdate(idx)     = idxIter;
     fD.deltaLL(idx)       = deltaLL(idxIter);
  end
  fD.Efy = Efy;   % what we did achieve in quality by the end
  lambdaHat = fD.lambdaTrace(:,idx);
  load(['/home/marcel/criticalityIterScaling/results/',runIdx, '.mat']) % load pars, fitoptions etc. to add to final file
  disp('Storing results up to (n-1)-ith step')     
     save(['results_',runIdx,'.mat'],...
           'fitoptions', 'pars', 'beta', 'runIdx', ...  % simulation setup
           'lambdaTrue', 'mfxTrain', ...                % desired outcome
           'lambdaHat', 'fD')                           % actual outcome
       
     j = size(diri,1);  
     idx = iteri(j); % idx may differ from j if some files are missing
     fname = diri(j,:); fname(fname==' ') = [];
      disp([ 'loading data ', num2str(j), '/', num2str(size(diri,1))] )
     load(['/home/marcel/criticalityIterScaling/results/',fname]);         % now load stuff and collect!
     fD.deltaLLs(:,idx)    = deltaLL;
     fD.deltas(:,idx)      = deltaIter;
     fD.lambdaTrace(:,idx) = lambdaIter;
     fD.EfyTrace(:,idx)    = Efy;
     fD.x0(:,idx)          = x0Iter;
     fD.idxUpdate(idx)     = idxIter;
     fD.deltaLL(idx)       = deltaLL(idxIter);
     
     fD.Efy = Efy;   % what we did achieve in quality by the end
     lambdaHat = fD.lambdaTrace(:,end);
     load(['/home/marcel/criticalityIterScaling/results/',runIdx, '.mat']) % load pars, fitoptions etc. to add to final file
     disp('Storing results')     
     save(['results_',runIdx,'.mat'],...
           'fitoptions', 'pars', 'beta', 'runIdx', ...  % simulation setup
           'lambdaTrue', 'mfxTrain', ...                % desired outcome
           'lambdaHat', 'fD')                           % actual outcome
     clear fitoptions pars beta runIdx lambdaTrue mfxTrain lambdaHat fD
     clear idx diri
end % end for i = runIDs

end