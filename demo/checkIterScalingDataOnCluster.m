
idxRun = [1:8];
d = 100; 
idxTouched = 1:d;

load('C:\Users\Loki\Desktop\tmp\timeStampsClusterRuns_bwVK.mat')


ifFigOverviewIndivResults = true; % make overview figure of 
 timePoint = 0.5;  % time point in days for which to illustrate results

ifFigurePerformanceTraces = true;

lgnd = {'800', '1600', '3200', '6400', '12800', ...
        '400 2^{(i/2000)}', '800 2^{(i/2000)}', '800 2^{(i/1000)}'};
    
nSmooth = 20;
nTimePoints = 1000;
timePoints = linspace( 1/24/60 , 3, nTimePoints);

%% Do necessary pre-computations

fDescr=[[1:d; (1:d)*nan],nchoosek(1:d,2)',[0:d; (0:d)*nan]]; 
 % key index for relating parameters lambda to expected values E[f(X)]
 
disp('Extracting time stamps')
times  = cell(length(idxRun),1);
fnames = cell(length(idxRun),1);
for i = idxRun
  tmp = dates(idxSess{i});
  tmp = tmp(2:end, :);   % remove sess*_***.mat, i.e. base file.
  times{i} = datenum(tmp(1:10:end, :));% thin out to match s*_res_small.mat
  tmp = names(idxSess{i});
  tmp = tmp(2:end, :);   % remove sess*_***.mat, i.e. base file.
  fnames{i} = tmp(1:10:end, :);        % thin out to match s*_res_small.mat
  times{i} = times{i} - times{i}(1,:);
end

clrs = hsv(length(idxRun));

% get time points
idxTimes = zeros(length(idxRun),nTimePoints);
for j = 1:nTimePoints
 for i = idxRun 
  [~, idxTimes(i,j)] = min(abs(times{i}-timePoints(j)));   
 end
end
idxTime = zeros(length(idxRun),1); % single chosen point for overview plots
for i = idxRun 
 [~, idxTime(i)] = min(abs(times{i}-timePoint));   
end

disp('Computing MSEs')

 idxK = 9:20;
 
% compute mean squared errors for expected values E[f(X)] and parameters 
MSEf = zeros(length(idxRun), nTimePoints);
MSEl = zeros(length(idxRun), nTimePoints);
MSEfh = zeros(length(idxRun), nTimePoints);
MSElh = zeros(length(idxRun), nTimePoints);
MSEfJ = zeros(length(idxRun), nTimePoints);
MSElJ = zeros(length(idxRun), nTimePoints);
MSEfV = zeros(length(idxRun), nTimePoints);
MSElV = zeros(length(idxRun), nTimePoints);
MSElVrw = zeros(length(idxRun), nTimePoints); % re-weighted by importance
for i = 1:length(idxRun)
 kRun = idxRun(i);
 load(['C:\Users\Loki\Desktop\tmp\s', num2str(kRun), '_res_small.mat'])
 MSEf(i, :) = mean( bsxfun(@minus, out.Efx, out.Efy(:,idxTimes(i,:))).^2 );
 MSEl(i, :) = mean( bsxfun(@minus, out.lambdaTrue, out.lambda(:,idxTimes(i,:))).^2 );
 MSEfh(i, :) = mean( bsxfun(@minus, out.Efx(1:d), out.Efy(1:d,idxTimes(i,:))).^2 );
 MSElh(i, :) = mean( bsxfun(@minus, out.lambdaTrue(1:d), out.lambda(1:d,idxTimes(i,:))).^2 );
 MSEfJ(i, :) = mean( bsxfun(@minus, out.Efx(d+1:d*(d+1)/2), out.Efy(d+1:d*(d+1)/2,idxTimes(i,:))).^2 );
 MSElJ(i, :) = mean( bsxfun(@minus, out.lambdaTrue(d+1:d*(d+1)/2), out.lambda(d+1:d*(d+1)/2,idxTimes(i,:))).^2 );
 
 MSEfV(i, :) = mean( bsxfun(@minus, out.Efx(end-d+idxK), out.Efy(end-d+idxK,idxTimes(i,:))).^2 );
 MSElV(i, :) = mean( bsxfun(@minus, out.lambdaTrue(end-d+idxK), out.lambda(end-d+idxK,idxTimes(i,:))).^2 );
 
 tmp = bsxfun(@minus, out.lambda(end-d+idxK,idxTimes(i,:)), out.lambdaTrue(end-d+idxK));
 covfk = diag(out.Efx(end-d+idxK)) - out.Efx(end-d+idxK) * out.Efx(end-d+idxK)';
% covfk = diag(ones(length(idxK),1));
 for j = 1:size(MSElVrw, 2)
  MSElVrw(i,j) = tmp(:,j)' * covfk * tmp(:,j);  
 end
 MSElVrw(i,:) = MSElVrw(i,:) / length(idxK); 
end


% get same stuff, but smoothed 
disp('Computing smoothed MSEs')
MSEfs = zeros(length(idxRun), nTimePoints-nSmooth);
MSEls = zeros(length(idxRun), nTimePoints-nSmooth);
MSEfhs = zeros(length(idxRun), nTimePoints-nSmooth);
MSElhs = zeros(length(idxRun), nTimePoints-nSmooth);
MSEfJs = zeros(length(idxRun), nTimePoints-nSmooth);
MSElJs = zeros(length(idxRun), nTimePoints-nSmooth);
MSEfVs = zeros(length(idxRun), nTimePoints-nSmooth);
MSElVs = zeros(length(idxRun), nTimePoints-nSmooth);
MSElVrws = zeros(length(idxRun), nTimePoints-nSmooth);

for i = 1:length(idxRun)
  disp(['i = ', num2str(i),'/',num2str(length(idxRun))])
 for j = 1:nTimePoints-nSmooth
  MSEfs(i,j) = mean(MSEf(i,j:j+nSmooth));
  MSEls(i,j) = mean(MSEl(i,j:j+nSmooth));
  MSEfhs(i,j) = mean(MSEfh(i,j:j+nSmooth));
  MSEfJs(i,j) = mean(MSEfJ(i,j:j+nSmooth));
  MSEfVs(i,j) = mean(MSEfV(i,j:j+nSmooth));
  MSElhs(i,j) = mean(MSElh(i,j:j+nSmooth));
  MSElJs(i,j) = mean(MSElJ(i,j:j+nSmooth));
  MSElVs(i,j) = mean(MSElV(i,j:j+nSmooth));
  MSElVrws(i,j) = mean(MSElVrw(i,j:j+nSmooth));
 end
end

%%

% Compute correlation matrices for chosen fixed time point
disp('Computing cij for fixed chosen time points')
mfxEval  = zeros(5151,length(idxRun)); % = E_hat[f(X)]
mfxTrain = zeros(5151,length(idxRun)); % = E_emp[f(X)]
Ceval   = cell(length(idxRun),1);
Ctrain  = cell(length(idxRun),1);
cLow = zeros(length(idxRun),1);
cHigh = zeros(length(idxRun),1);

for k = 1:length(idxRun);
    
kRun = idxRun(k);
load(['C:\Users\Loki\Desktop\tmp\s', num2str(kRun), '_res_small.mat'])
mfxEval(:,k) = out.Efy(:,idxTime(kRun)); % pick result at chosen time point
 %load(['C:\Users\Loki\Desktop\tmp\EfxHat']) % alternatively, pick                                             
 %mfxEval = EfxHat(:,kRun);                  % well-sampled final results
mfxTrain(:,k) = out.Efx;
clear out

Ceval{k} = zeros(d,d);
Ctrain{k} = zeros(d,d);
ij = zeros(d,d);
for i = 1:d      % off-diagonal 
 for j = (i+1):d % entries
   ij(i,j) = find((fDescr(1,:)==i)&(fDescr(2,:)==j),1,'first');
   Ceval{k}(i,j)  = mfxEval( ij(i,j),k ) - mfxEval(i,k)*mfxEval(j,k) ;
   Ctrain{k}(i,j) = mfxTrain(ij(i,j),k)  - mfxTrain(i,k)*mfxTrain(j,k) ; 
 end
end
for i = 1:d % diagonal entries
   Ceval{k}(i,i) = mfxEval(i,k) * (1- mfxEval(i,k));
   Ctrain{k}(i,i) = mfxTrain(i,k) * (1- mfxTrain(i,k));
end
for i = 1:d % make covariance matrices to be correlation matrices
 for j = i:d
   Ceval{k}(i,j)  = Ceval{k}(i,j)  / sqrt( Ceval{k}(i,i) ) / sqrt( Ceval{k}(j,j) );
   Ctrain{k}(i,j) = Ctrain{k}(i,j) / sqrt( Ctrain{k}(i,i)) / sqrt(Ctrain{k}(j,j) ); 
 end
end
Ceval{k} = Ceval{k} + Ceval{k}' - diag(diag(Ceval{k}));     % make 
Ctrain{k} = Ctrain{k} + Ctrain{k}' - diag(diag(Ctrain{k})); % symmetric

cLow(k) = min([Ctrain{k}(:);Ceval{k}(:)]); % info for proper visualization
cHigh(k) = max([Ctrain{k}(Ctrain{k}<0.999);Ceval{k}(Ceval{k}<0.999)]);

end

% compute time traces for covariances / correlations
disp('Computing MSE in cij for all chosen time points')
MSEcov = zeros(length(idxRun), nTimePoints);
MSEcor = zeros(length(idxRun), nTimePoints);
covij = cell(length(idxRun),1);
corij = cell(length(idxRun),1);
covTrue = zeros(d*(d-1)/2, length(idxRun));
corTrue = zeros(d*(d-1)/2, length(idxRun));
for i = 1:length(idxRun)
disp(['i = ', num2str(i),'/',num2str(length(idxRun))])
kRun = idxRun(i);
load(['C:\Users\Loki\Desktop\tmp\s', num2str(kRun), '_res_small.mat'])
  covij{i} = out.Efy((d+1 : d*(d+1)/2),idxTimes(i,1:nTimePoints)) ...
        - (out.Efy(fDescr(1, (d+1 : d*(d+1)/2)),idxTimes(i,1:nTimePoints)) ... 
        .* out.Efy(fDescr(2, (d+1 : d*(d+1)/2)),idxTimes(i,1:nTimePoints)));
  corij{i} = covij{i} ...
    ./ sqrt( out.Efy(fDescr(1, (d+1 : d*(d+1)/2)),idxTimes(i,1:nTimePoints))  ...
         .* ( 1- out.Efy(fDescr(1, (d+1 : d*(d+1)/2)),idxTimes(i,1:nTimePoints)) )) ...
    ./ sqrt( out.Efy(fDescr(1, (d+1 : d*(d+1)/2)),idxTimes(i,1:nTimePoints)) ...
          .* ( 1- out.Efy(fDescr(1, (d+1 : d*(d+1)/2)),idxTimes(i,1:nTimePoints)) ));
  covTrue(:,i) = out.Efx((d+1 : d*(d+1)/2)) ...
            - (out.Efx(fDescr(1,(d+1 : d*(d+1)/2))) .* out.Efx(fDescr(2,(d+1 : d*(d+1)/2))));
  varTrue = out.Efx(1:d) .* ( 1- out.Efx(1:d));
  corTrue(:,i) = covTrue(:,i) ./ sqrt(varTrue(fDescr(1, (d+1 : d*(d+1)/2)))) ./ sqrt(varTrue(fDescr(2, (d+1 : d*(d+1)/2))));
  MSEcov(i, :) = mean( bsxfun(@minus, covTrue(:,i), covij{i}).^2 );
  MSEcor(i, :) = mean( bsxfun(@minus, corTrue(:,i), corij{i}).^2 );
  covij{i} = []; % have to save 
  corij{i} = []; % on memory
end
clear covij corij ans compile_c i j k kRun tmp varTrue

allVKs = zeros(length(idxRun), (d+1), nTimePoints);
for i = 1:length(idxRun)
  kRun = idxRun(i);
  load(['C:\Users\Loki\Desktop\tmp\s', num2str(kRun), '_res_small.mat'])
  allVKs(i, :, :) = out.lambda(end-d:end,idxTimes(i,:));
end


















%%
if ifFigOverviewIndivResults
%--------------------------------------------------------------------------
% Create first overview figure with
% a) E[x_i]            true vs model fit
% b) E[x_i  * x_j]     true vs model fit
% c) E[K == k]         
% d) Correlations c_ij true vs model fit
% e) data correlation matrix
% f) fit correlation matrix

for kRun = idxRun
    
figure; 
subplot(3,2,1)
line([0,1],[0,1], 'color', 'k', 'linestyle', '--')
hold on
plot(mfxTrain(1:d,k), mfxEval(1:d,k), 'k.', 'markerSize', 5) 
xlabel('Data'); ylabel('Model draw')
title('E[x_i]')
subplot(3,2,2)
line([0,1],[0,1], 'color', 'k', 'linestyle', '--')
hold on
plot(mfxTrain(d+1:end-d-1,k), mfxEval(d+1:end-d-1,k), 'k.','markerSize', 5) 
%idxTouchedJ = idxTouched; 
%idxTouchedJ(idxTouchedJ<d) = [];
%idxTouchedJ(idxTouchedJ>d*(d+1)/2) = [];
%plot(mfxTrain(idxTouchedJ), mfxEval(idxTouchedJ), 'r.', 'markerSize', 5)
xlabel('Data'); ylabel('Model draw')
axis([0,1,0,1])
title('E[x_i x_j]')

subplot(3,2,3)
line([cLow(k),cHigh(k)],[cLow(k),cHigh(k)], 'color','k','linestyle', '--')
hold on
plot(Ctrain{k}(:), Ceval{k}(:), 'k.', 'markerSize', 5) 
%idxTouchedCorr = false(d,d);
%for i = 1:d
%    for j = 1:d
%        idxTouchedCorr(i,j) = (ismember(ij(i,j), idxTouched)==1);
%    end
%end
%plot(Ctrain{k}(idxTouchedCorr), Ceval{k}(idxTouchedCorr), 'r.', 'markerSize', 5)
%for i = 1:d
% if ~ismember(i,idxTouched(1:d))
%  idxTouchedCorr(i,:) = false;
%  idxTouchedCorr(:,i) = false;
% end
%end
%plot(Ctrain{k}(idxTouchedCorr), Ceval{k}(idxTouchedCorr), 'g.', 'markerSize', 15)
xlabel('Data'); ylabel('Model draw')
axis([cLow(k),cHigh(k),cLow(k),cHigh(k)])
title('corr(x_i,x_j)')


subplot(3,2,5)
imagesc(Ctrain{k}-diag(diag(Ctrain{k}))); 
set(gca, 'CLim', [cLow(k), cHigh(k)]);
title('correlation matrix data')
subplot(3,2,6)
imagesc(Ceval{k}-diag(diag(Ceval{k}))); 
set(gca, 'CLim', [cLow(k), cHigh(k)]);
title('correlation matrix model draw')
subplot(3,2,4)
line([0,1],[0,1], 'color', 'k', 'linestyle', '--')
hold on
plot(mfxTrain(end-d:end,k),mfxEval(end-d:end,k), 'k.', 'markerSize', 5)
VHigh = max([mfxTrain(end-d:end,k);mfxEval(end-d:end,k)]);
xlabel('Data'); ylabel('Model draw')
axis([0,VHigh,0,VHigh])
title('E[ \Sigma_{i=1}^n x_i = K ]')

ttl = ['run number ', num2str(kRun)];
annotation('textbox', [0 0.9 1 0.1], ...
    'String', ttl, ... 
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')   

end
end % end if ifFigOverviewIndivResults

%%
if ifFigurePerformanceTraces

figure; 
subplot(1,2,1)
for i = 1:length(idxRun)
 plot(timePoints, log(MSEf(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('Mean squared error of target expected values of features E[f(X)]')

subplot(1,2,2)
for i = 1:length(idxRun)
 plot(timePoints, (MSEl(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
legend(lgnd(idxRun));
legend boxoff
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('Mean squared error of target parameter vectors \lambda_{True}')

%%
figure; 
subplot(2,3,1)
for i = 1:length(idxRun)
 plot(timePoints, log(MSEfh(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for features E[x_i] (i.e. for h)')

subplot(2,3,2)
for i = 1:length(idxRun)
 plot(timePoints, log(MSEcov(i,:)),'color',clrs(i,:),'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for covariances cov(x_i, x_j) (corresponding to J)')

subplot(2,3,3)
for i = 1:length(idxRun)
 plot(timePoints, log(MSEfV(i,:)),'color',clrs(i,:),'linewidth',2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for features E[K = k] (i.e. for V(K))')

subplot(2,3,4)
for i = 1:length(idxRun)
 plot(timePoints, (MSElh(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters h')

subplot(2,3,5)
for i = 1:length(idxRun)
 plot(timePoints, (MSElJ(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters J')

subplot(2,3,6)
for i = 1:length(idxRun)
 plot(timePoints, (MSElV(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters V(K)')

%%
figure; 
subplot(2,3,1)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), log(MSEfhs(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for features E[x_i] (i.e. for h)')

subplot(2,3,2)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), log(MSEfJs(i,:)),'color',clrs(i,:),'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for features E[x_i x_j] (i.e. for J)')

subplot(2,3,3)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), log(MSEfVs(i,:)),'color',clrs(i,:),'linewidth',2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for features E[K = k] (i.e. for V(K))')

subplot(2,3,4)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), (MSElhs(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters h')

subplot(2,3,5)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), (MSElJs(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters J')

subplot(2,3,6)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), (MSElVs(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters V(K)')

%% MSE for covariances 
figure; 
subplot(2,2,1)
for i = 1:length(idxRun)
 plot(timePoints, log(MSEcov(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for covariances')

subplot(2,2,2)
for i = 1:length(idxRun)
 plot(timePoints, log(MSEcor(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('log(MSE)')
title('log(MSE) for correlations')

subplot(2,2,3)
for i = 1:length(idxRun)
 plot(timePoints, (MSElVrw(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for re-weighted parameters \lambda for terms V(K)')

subplot(2,2,4)
for i = 1:length(idxRun)
 plot(timePoints(1:end-nSmooth), (MSElVrws(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('smoothed MSE for re-weighted parameters \lambda for terms V(K)')

%%

figure;
subplot(2,3,1)
plot(0:50, out.Efx(end-d+(0:50)), 'ko-', 'linewidth', 2);
box off
xlabel('K')
ylabel('P(K)')
title('Empirical distribution P(K=k)')

subplot(2,3,2)
for i = 1:length(idxRun)
 plot(timePoints, (MSElV(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for parameters V(K)')

subplot(2,3,5)
for i = 1:length(idxRun)
 plot(timePoints, (MSElVrw(i,:)), 'color', clrs(i,:), 'linewidth', 2);
 hold on
end
box off
set(gca, 'TickDir', 'out')
xlabel('computation time (in days)')
ylabel('MSE')
title('MSE for re-w. parameters \lambda for terms V(K), 8 < K < 20')

subplot(2,3,6)
plot(out.lambdaTrue(end-d+(0:50)), 'ko-', 'linewidth', 2)
hold on
plot(out.fitoptions.lambda0(end-d+(0:50)), 'bo--', 'linewidth', 2)
legend('True parameters', 'Initial guess')
title('True \lambda and fitting initialization \lambda_0')
for i = 1:8
  plot(squeeze(allVKs(i,1:51,end)), '.-', 'color', clrs(i,:))
end
box off
set(gca, 'TickDir', 'out')
xlabel('K')
ylabel('V(K)')
title('Parameters \lambda, true, initial and final fits')












%% TEST SOME STUFF
for i = 1:length(idxRun)
 kRun = idxRun(i); 
 load(['C:\Users\Loki\Desktop\tmp\s', num2str(kRun), '_res_small.mat'])
 idxK = 0:25;
 figure; 
 subplot(3,1,1); 
 imagesc(diff(out.lambda(end-100+idxK,idxTimes(kRun,:)),1,2)); colorbar; 
 for i = 1:3
  [~, xlabelTimePoints(i)] = min(abs(timePoints-i));
 end
% set(gca, 'XTick', idxTimes(kRun, xlabelTimePoints));
% set(gca, 'XTickLabel', 1:3);
 xlabel('K')
 xlabel('time in days')
 title('changes in entries of \lambda corresponding to V(K)')
 subplot(3,1,2); 
 imagesc(bsxfun(@minus, out.lambda(end-100+idxK,idxTimes(kRun,:)), out.lambdaTrue(end-100+idxK))); colorbar;
% set(gca, 'XTick', idxTimes(kRun, xlabelTimePoints));
% set(gca, 'XTickLabel', 1:3);
 xlabel('K')
 xlabel('time in days')
 title('absolute E[f(X)] for features f corresponding to V(K)')
 subplot(3,1,3); 
% imagesc( diff(out.Efy(end-100+idxK,idxTimes(kRun,200:end)),1,2) ); colorbar;
 imagesc(out.Efy(end-100+idxK,idxTimes(kRun,:))); colorbar;
% set(gca, 'XTick', idxTimes(kRun, xlabelTimePoints));
% set(gca, 'XTickLabel', 1:3);
 xlabel('K')
 xlabel('time in days')
 title('changes in E[f(X)] for features corresponding to V(K)')
 
 pause; 
end


%%

figure; 
for k = 2:d+1
 subplot(4,5,k-1)
 for i = 1:length(idxRun)
   kRun = idxRun(i);
   plot(timePoints(1:250), squeeze(allVKs(i,k,1:250)), 'color', clrs(i,:));
   hold on
 end
   plot(0.8, out.lambdaTrue(end-d+(k-1)), '*', 'color', 'g');
   plot(0.8, out.fitoptions.lambda0(end-d+(k-1)), '*', 'color', 'r');
   %title(['Real (g*), initial (r*) and progessively estimated (traces) V(K) , K =', num2str(k)-1])
   if mod(k-1,5)==1
    ylabel('V(K)')
   end
   if k-1 > 15
    xlabel('time in days')
    set(gca, 'XTick', [0, 0.5])
   else 
    set(gca, 'XTick', []);
   end
   box off
   axis([0, 0.85, 0, 1])
   axis autoy
 pause;
 hold off
end
%%
figure;
i = 5;
for j = 1:50:1000, 
 plot(cTrue(:,i), covij{i}(:, j), '.'); 
 line([-0.02, 0.1], [-0.02, 0.1], 'color', 'k'); 
 title(['step #', num2str(10*idxTimes(i,j)), ', time is ', num2str(timePoints(j)) ' days']), 
 pause; 
end
end % end if ifFigurePerformanceTraces