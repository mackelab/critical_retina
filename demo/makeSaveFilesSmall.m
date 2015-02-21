function makeSaveFilesSmall(fname, n, thinning)
n = 10*n;

if nargin < 3 || isempty(thinning)
    thinning = 10;
end

load(['/home/marcel/criticalityIterScaling/results/', fname, 'n', num2str(n), 'asfinal.mat'])

for i = 1:length(fD)
 if ~isempty(fD{i})
  fD{i}.deltaLLs    = fD{i}.deltaLLs(:,1:thinning:end);
  fD{i}.deltas      = fD{i}.deltas(:,1:thinning:end);
  fD{i}.lambdaTrace = fD{i}.lambdaTrace(:,1:thinning:end);
  fD{i}.EfyTrace    = fD{i}.EfyTrace(:,1:thinning:end);
  fD{i}.x0          = fD{i}.x0(:,1:thinning:end);
  fD{i}.idxUpdate   = fD{i}.idxUpdate(:,1:thinning:end);
  fD{i}.dletaLLs    = fD{i}.deltaLL(:,1:thinning:end);
  fD{i}.MSEperc     = fD{i}.MSEperc(:,1:thinning:end);
  fD{i}.MSE         = fD{i}.MSE(:,1:thinning:end); 
 end
end
save(['/home/marcel/criticalityIterScaling/results/', fname, 'n', num2str(n), 'assmall.mat'], 'lambdaHat','fD', 'thinning')

end