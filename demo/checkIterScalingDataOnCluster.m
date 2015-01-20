
idxRun = [1:8];
d = 100; 
idxTouched = 1:d;

% Create first overview figure with
% a) E[x_i]            true vs model fit
% b) E[x_i  * x_j]     true vs model fit
% c) E[K == k]         
% d) Correlations c_ij true vs model fit
% e) data correlation matrix
% f) fit correlation matrix
for kRun = idxRun
load(['C:\Users\Loki\Desktop\tmp\s', num2str(kRun), '_res_small.mat'])
mfxEval = out.Efy(:,end);
mfxTrain = out.Efx;
clear out

figure; 
subplot(3,2,1)
line([0,1],[0,1], 'color', 'k', 'linestyle', '--')
hold on
plot(mfxTrain(1:d), mfxEval(1:d), 'k.', 'markerSize', 5) 
xlabel('Data'); ylabel('Model draw')
title('E[x_i]')

fDescr=[[1:d; (1:d)*nan],nchoosek(1:d,2)',[0:d; (0:d)*nan]];
Ceval = zeros(d,d);
Ctrain = zeros(d,d);
ij = zeros(d,d);
for i = 1:d
 for j = (i+1):d
   ij(i,j) = find((fDescr(1,:)==i)&(fDescr(2,:)==j),1,'first');
   Ceval(i,j)  = mfxEval( ij(i,j) ) - mfxEval(i)*mfxEval(j) ;
   Ctrain(i,j) = mfxTrain(ij(i,j))  - mfxTrain(i)*mfxTrain(j) ; 
 end
end
for i = 1:d
   Ceval(i,i) = mfxEval(i) * (1- mfxEval(i));
   Ctrain(i,i) = mfxTrain(i) * (1- mfxTrain(i));
end
for i = 1:d
 for j = i:d
   Ceval(i,j)  = Ceval(i,j)  / sqrt( Ceval(i,i) ) / sqrt( Ceval(j,j) );
   Ctrain(i,j) = Ctrain(i,j) / sqrt( Ctrain(i,i)) / sqrt(Ctrain(j,j) ); 
 end
end
Ceval = Ceval + Ceval' - diag(diag(Ceval));
Ctrain = Ctrain + Ctrain' - diag(diag(Ctrain));

cLow = min([Ctrain(:);Ceval(:)]); cHigh = max([Ctrain(Ctrain<0.999);Ceval(Ceval<0.999)]);

subplot(3,2,2)
line([0,1],[0,1], 'color', 'k', 'linestyle', '--')
hold on
plot(mfxTrain(d+1:end-d-1), mfxEval(d+1:end-d-1), 'k.', 'markerSize', 5) 
idxTouchedJ = idxTouched; 
idxTouchedJ(idxTouchedJ<d) = [];
idxTouchedJ(idxTouchedJ>d*(d+1)/2) = [];
%plot(mfxTrain(idxTouchedJ), mfxEval(idxTouchedJ), 'r.', 'markerSize', 5)
xlabel('Data'); ylabel('Model draw')
axis([0,1,0,1])
title('E[x_i x_j]')

subplot(3,2,3)
line([cLow,cHigh],[cLow,cHigh], 'color', 'k', 'linestyle', '--')
hold on
plot(Ctrain(:), Ceval(:), 'k.', 'markerSize', 5) 
idxTouchedCorr = false(d,d);
for i = 1:d
    for j = 1:d
        idxTouchedCorr(i,j) = (ismember(ij(i,j), idxTouched)==1);
    end
end
plot(Ctrain(idxTouchedCorr), Ceval(idxTouchedCorr), 'r.', 'markerSize', 5)
for i = 1:d
 if ~ismember(i,idxTouched(1:d))
  idxTouchedCorr(i,:) = false;
  idxTouchedCorr(:,i) = false;
 end
end
plot(Ctrain(idxTouchedCorr), Ceval(idxTouchedCorr), 'g.', 'markerSize', 15)
xlabel('Data'); ylabel('Model draw')
axis([cLow,cHigh,cLow,cHigh])
title('corr(x_i,x_j)')


subplot(3,2,5)
imagesc(Ctrain-diag(diag(Ctrain))); 
set(gca, 'CLim', [cLow, cHigh]);
title('correlation matrix data')
subplot(3,2,6)
imagesc(Ceval-diag(diag(Ceval))); 
set(gca, 'CLim', [cLow, cHigh]);
title('correlation matrix model draw')
subplot(3,2,4)
line([0,1],[0,1], 'color', 'k', 'linestyle', '--')
hold on
plot(mfxTrain(end-d:end),mfxEval(end-d:end), 'k.', 'markerSize', 5)
VHigh = max([mfxTrain(end-d:end);mfxEval(end-d:end)]);
xlabel('Data'); ylabel('Model draw')
axis([0,VHigh,0,VHigh])
title('E[ \Sigma_{i=1}^n x_i = K ]')

ttl = ['run number ', num2str(kRun)];
annotation('textbox', [0 0.9 1 0.1], ...
    'String', ttl, ... 
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')   

end