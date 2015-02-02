function getEfxHat(fname, nSamples, burnIn)

if nargin < 3 || isempty(burnIn)
    burnIn = 16000;
end
if nargin < 2 || isempty(nSamples)
    nSamples = 32000;
end

folder = '/home/marcel/criticalityIterScaling/results/';
if exist([folder,fname], 'file')
  load([folder, fname]);
   % gives lambdaHat, fD (fit diagnostics), fitoptions, beta, n , idxRepet
else % assemble data from currently most recent iteration step
  warning('assembling data from unfinished runs not yet implemented')  
end

EfxHat = zeros(n*(n+3)/2 +1, max(idxRepet)); 
EE     = zeros(  2,          max(idxRepet));
varE   = zeros(  1,          max(idxRepet));
for i = 1:length(idxRepet)
  j = idxRepet(i);
  if ~ismember(j, idxRepet)
     error('Results for the desired data set were not stored in this file')
  end
  [EfxHat(:,j), EE(:,j), ~] = maxEnt_gibbs_pair_C(nSamples, burnIn, lambdaHat(:,j), n, 'cluster');
  varE(j) = EE(2,j) - EE(1,j)^2;
end
save([folder, fname], 'lambdaHat', 'fD', 'fitoptions', 'beta', 'n', ...
                      'idxRepet', 'EfxHat', 'EE', 'varE');

