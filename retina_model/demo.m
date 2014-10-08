
d = [100, 150]; 
n = 30; 

mode.tiling = 'random';
mode.RF = 'center-surround DoG';

pars.thres = 0.00001;

Sigma = { [15,   0 ;      % template for central part of DoG filters
            0,  15 ], ... % (peak of Mexican heat)

          [20,   0 ;      % template for outer part of DoG filters
            0,  20 ]};    % (surround of Mexican hat)
        
SigmaDFs = [1000, 1000];  % degrees of freeom for DoG component covariances

ON_OFF = 2*randi(2, [n,1]) - 3; % per cell whether it's an ON or an OFF RGC

hight = [0.9, 1]; % template for hights of central and outer DoG components

hightSTDs = [0.01, 0.01]; % standard deviates for hights of DoG components 

idxON  = find(ON_OFF > 0);
idxOFF = find(ON_OFF < 0);
pars.hight = zeros(n,2);
pars.hight(:,1) = abs(normrnd(hight(1),hightSTDs(1)^2,[n,1])) .*ON_OFF;
pars.hight(:,2) = abs(normrnd(hight(2),hightSTDs(2)^2,[n,1])) .*ON_OFF;
for i = 1:n
  pars.Sigma{i,1} = wishrnd(Sigma{1}, SigmaDFs(1))/SigmaDFs(1);
  pars.Sigma{i,2} = wishrnd(Sigma{2}, SigmaDFs(2))/SigmaDFs(2);
end

%%-------------------------------------------------------------------------

[W, RGCcen] = genFilters(d,n,mode,pars); 

%%-------------------------------------------------------------------------

figure(1); 
subplot(1,2,1)
for i = 1:length(idxON); % ON cells 
 xy = reshape(W(idxON(i),:), d(1), d(2)); 
 xy(xy==0) = NaN;                  % don't plot these pixels, just clutters
 surf(xy, ones(d)*i/length(idxON)), hold on, 
end;
hold off
title('ON RGCs')

subplot(1,2,2)
for i = 1:length(idxOFF); % OFF cells
 xy = reshape(W(idxOFF(i),:), d(1), d(2)); 
 xy(xy==0) = NaN;                  % don't plot these pixels, just clutters
 surf(xy, ones(d)*i/length(idxOFF)), hold on, 
end;
hold off
title('OFF RGCs')

figure(2); 
subplot(1,2,1)
xy = reshape(sum(W(idxON,:),1), d(1), d(2));
imagesc(xy),  
title('ON RGCs')
subplot(1,2,2)
xy = reshape(sum(W(idxOFF,:),1), d(1), d(2));
imagesc(xy),  
title('OFF RGCs')
