%% produces fig 1 for the journal version of 'criticality and correlations' 
% figure '2' is the overview figure summarising flat models, how we fit 
% them to data and that we still obtain signatures of 'criticality' when
% applying them to our data

clear all

ifAddGaussianFits = false;

splitFigs = false;

if ~splitFigs
  figure2 = figure('Tag', 'fig2', 'units','normalized', ...
                   'position',[0,0,0.99,0.99]);
end

addpath(genpath('../code/'))

axesThickness  = 1.1; % overall thickness of axes. Doesn't seem to do much.

fontName = 'Arial';    fontWeight     = 'normal';
fontSize       = 1 * 10;   fontSizeTitle  = 1 * 16;   
fontSizeXlabel = 1 * 10;   fontSizeYlabel = 1 * 11;
fontSizeText   = 1 * 10;   fontSizeLegend = 1 * 11;

load('../data/RGC_sim_nat.mat')        % loads raw simulated RGC data

 n = size(output,3);         % get input spike raster format
 Nc = size(output,1);        % n = number of cells, Nc = trial length
 nTrials = size(output,2);   % nTrials = number of trials
 output = double(output);      % pyhton code produces data of type int64
 tmp = zeros(n, Nc, nTrials);
 for i = 1:nTrials
   tmp(:, :, i) = squeeze(output(:,i,:))';
 end
clear output;     % could pick a new name, but I got kind of to this one 
output = struct;  %
output.spikes = zeros(n, Nc*nTrials);
for i = 1:n
   tmptmp = squeeze(tmp(i,:,:));
   output.spikes(i, :) = tmptmp(:);
end
output.spkCorrs = corr(output.spikes');
clear i tmp tmptmp

maxN = 120;
maxi = maxN / 20;
clrs = [254,153,41;
236,112,20;
204, 76,2;      % colors for inidivual traces.
153,52,4;
102,37,6;
0,0,0]/255;

%% a) different model fits to simulated data
if splitFigs
  figure(21)
  subplot(1,2,1)
else
  subplot(21,34,vec(bsxfun(@plus,  (1:8)', (0:9)*34)))
end
 n = size(output.spikes,1);
 pcount = sum(output.spikes,1)';
 pcount = histc(pcount,0:n);
 pcount = pcount/sum(pcount);
 
  mu1 =  (0:n) * pcount;
  mu2 = (0:n).^2 * pcount ;
  Z = ( (n) * (mu2/mu1 - mu1 -1)) + mu1;
  a = (n * mu1 - mu2) / Z;
  b = (n - mu1) * (n - mu2/mu1) / Z; 
  lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';  
  pcount_betabin = exp(lognchoosek + betaln(a + (0:n), ...
                        n + b - (0:n))' - betaln(a, b)); 

 [mu,rho,~]=meanvar_count_2_meancorr(mu1,mu2 - mu1^2,n); 
 [hFI,JFI,pcount_flatIsing,~]=fit_flat_ising_model(mu,rho,n);
 
 pcount_binomial = binopdf(0:n, n, mean(output.spikes(:)));

semilogy(0:n,pcount, 's-', 'color', 'k', 'linewidth', 1.5,'markerSize',1.5)
hold on
semilogy(0,pcount_betabin(1),   'color', clrs(2,:), ...
         'linestyle', '-', 'linewidth', 1.5);
semilogy(0:n,pcount_binomial,  'color', [72, 60, 50]/255, ...
         'linestyle', ':', 'linewidth', 1.5);
semilogy(0:n,pcount_flatIsing, 'color', [145, 163, 176]/255, ...
         'linestyle', '--', 'linewidth', 1.5);
semilogy(0:n,pcount, 's', 'color', 'k', 'linewidth', 0.5, 'markerSize',1.5) 
semilogy(0:n,pcount_betabin,   'color', clrs(2,:), 'linestyle', '-',  ...
         'linewidth', 2.5);
if ifAddGaussianFits 
 pars0 = [0;0];
 fitoptions = [];
 fitoptions.DerivativeCheck = 'on';
 fun = @(pars) computeGradient2oMaxEnt(pars, n, mu1, mu2);
 [parsOut, fval] = minFunc(fun, pars0, fitoptions);
 [f,g,Ex,Ex2] = computeGradient2oMaxEnt(parsOut, n, mu1, mu2);
 disp('desired and true E(K) and E(K^2)')
 [mu1, Ex; mu2, Ex2],  pause; 
 pcount_gauss = exp(parsOut(1) * (0:n) + parsOut(2) * (0:n).^2);
 pcount_gauss = pcount_gauss/sum(pcount_gauss);
 semilogy(0:n,pcount_gauss,   'color', 'r', 'linestyle', '--', ...
          'linewidth', 2.5);
end
hold off

legend('data', 'beta-binomial', 'binomial', 'flat Ising')
set(gca, 'FontSize', fontSize)  
xlabel('population spike count', 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight )
axis([0, 140, 1.1 * 10^(-5), 1]); 
box off, set(gca, 'TickDir' ,'out')

if splitFigs
  figure(21)
  subplot(1,2,2)
else
  subplot(21,34,vec(bsxfun(@plus,  (1:8)', (12:20)*34)))
end

idxRep = 1; 

 for i =1:maxi
   semilogy(0,0, '-', 'color', clrs(i,:), 'linewidth', 2.5)
   hold on
   lgnd{i} = ['n = ', num2str(20*i)];
 end
 legend(lgnd);
 for i =1:maxi
     
   pcount = zeros(20*i+1,1);              
   for t = 1:idxRep
    idxN = randsample(n,20*i);    % pick subset of Ns(j) <= n neurons
    histCounts = full(sum(output.spikes(idxN,:),1));          
    pcount = pcount + histc(histCounts,0:20*i)';   
   end     
   pcount = pcount/sum(pcount);
    
  semilogy(0:length(pcount)-1,pcount, 's-', 'color', clrs(i,:), ...
           'linewidth', 1.5, 'markerSize', 1.5)
  hold on
  semilogy(0:length(pcount)-1,pcount, 's', 'color', clrs(i,:), ...
           'linewidth', 0.75, 'markerSize', 1.5)
  mu1 =  (0:20*i) * pcount;
  mu2 = (0:20*i).^2 * pcount ;
  Z = ( (20*i) * (mu2/mu1 - mu1 -1)) + mu1;
  a = (20*i * mu1 - mu2) / Z;
  b = (20*i - mu1) * (20*i - mu2/mu1) / Z; 
  lognchoosek = (gammaln(20*i+1) - gammaln((0:20*i)+1) - ...
                 gammaln(20*i+1-(0:20*i)))';  
  logpcount = lognchoosek + betaln(a + (0:20*i), 20*i + b - (0:20*i))' ...
              - betaln(a, b); 
  pcount = exp(logpcount);

  semilogy(0:length(pcount)-1, pcount, '-', 'color', clrs(i,:), ...
           'linewidth', 2.5)
  if ifAddGaussianFits
   pars0 = [0;0];
   fitoptions = [];
   fitoptions.DerivativeCheck = 'on';
   f = @(pars) computeGradient2oMaxEnt(pars, 20*i, mu1, mu2);
   [parsOut, fval] = minFunc(f, pars0, fitoptions);
   [f,g,Ex,Ex2] = computeGradient2oMaxEnt(parsOut, 20*i, mu1, mu2);
   disp('desired and true E(K) and E(K^2)')
   [mu1, Ex; mu2, Ex2],  pause; 
   logpcount_gauss = (parsOut(1) * (0:20*i) + parsOut(2) * (0:20*i).^2);
   pcount_gauss = exp(logpcount_gauss - max(logpcount_gauss));
   pcount_gauss = pcount_gauss/sum(pcount_gauss);
   semilogy(0:length(pcount)-1, pcount_gauss, ':', 'color', clrs(i,:), ...
           'linewidth', 2.5)
  end  
 end
hold off

set(gca, 'FontSize', fontSize)  
xlabel('population spike count', 'FontName', fontName, 'FontSize', ...
       fontSizeXlabel, 'FontWeight', fontWeight )
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight )
axis([0, 51, 1.1*10^(-4), 1]); 
box off, set(gca, 'TickDir' ,'out')


%% b) different fits to other people's data
if splitFigs
  figure(22)
else
  subplot(21,34,vec(bsxfun(@plus,  (11:18)', (12:20)*34)))
end
for i = 1:2:5
  semilogy(-1, 0, 's-', 'color', clrs(i,:), 'linewidth',1.5,'markerSize',3)
  hold on;
end
legend('Okun et al. 2012', 'Tkacik et al. 2013', 'Tkacik et al. 2012')
% start with Okun data (n=96. Watch out: not 96 neurons, but 96 tetrodes!)
n = 96;
data_1=load('../data/other_studies/figure_data/Okun_Counts.txt');
ks = round(data_1(:,1));
pcount = data_1(:,2); pcount = pcount/sum(pcount);
mu1 =  (ks') * pcount;               mu2 = (ks').^2 * pcount;
Z = ( n * (mu2/mu1 - mu1 -1)) + mu1;
a = (n * mu1 - mu2) / Z;             b = (n - mu1) * (n - mu2/mu1) / Z;
lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';
logpcount = lognchoosek + betaln(a + (0:n), n + b - (0:n))' - betaln(a, b);

[a2,b2] = fitBetabinML(pcount, n, 100000);
logpcount2 = lognchoosek + betaln(a2 + (0:n), n + b2 - (0:n))' ...
            - betaln(a2, b2);
semilogy(data_1(:,1), pcount, 's-', 'color', clrs(1,:), ...
         'linewidth', 1.5, 'markerSize', 1.5)
hold on
semilogy(data_1(:,1), pcount, 's', 'color', clrs(1,:), ...
         'linewidth', 0.75, 'markerSize', 3)
semilogy(0:n, exp(logpcount), '-', 'color', clrs(1,:), 'linewidth', 2.5)

if ifAddGaussianFits
 pars0 = [0;0];
 fitoptions = [];
 fitoptions.DerivativeCheck = 'on';
 f = @(pars) computeGradient2oMaxEnt(pars, n, mu1, mu2);
 [parsOut,~] = minFunc(fun, pars0, fitoptions);
  [f,g,Ex,Ex2] = computeGradient2oMaxEnt(parsOut,n, mu1, mu2);
  disp('desired and true E(K) and E(K^2)')
  [mu1, Ex; mu2, Ex2],  pause; 
  pcount_gauss = exp(parsOut(1) * (0:n) + parsOut(2) * (0:n).^2);
  pcount_gauss = pcount_gauss/sum(pcount_gauss);
 semilogy(0:n, pcount_gauss, ':', 'color', clrs(1,:), 'linewidth', 2.5)
end
clear data_1 data_2 mu1 mu2 a b Z lognchoosek logpcount ks pcount

% next add data from Tkacik et al. 2014
n = 120;
data_1=load('../data/other_studies/figure_data/Tkacik_2014_data.txt');
ks = round(data_1(:,1));
pcount = data_1(:,2); pcount = pcount/sum(pcount);
mu1 =  (ks') * pcount;               mu2 = (ks').^2 * pcount;
Z = ( n * (mu2/mu1 - mu1 -1)) + mu1;
a = (n * mu1 - mu2) / Z;             b = (n - mu1) * (n - mu2/mu1) / Z;
lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';
logpcount = lognchoosek + betaln(a + (0:n), n + b - (0:n))' - betaln(a, b);

[a2,b2] = fitBetabinML(pcount, n, 100000);
logpcount2 = lognchoosek + betaln(a2 + (0:n), n + b2 - (0:n))' ...
             - betaln(a2, b2);
semilogy(data_1(:,1), pcount, 's-', 'color', clrs(3,:), ...
         'linewidth', 1.5, 'markerSize', 1.5)
hold on
semilogy(data_1(:,1), pcount, 's', 'color', clrs(3,:), ...
         'linewidth', 0.75, 'markerSize', 3)
semilogy(0:n, exp(logpcount), '-', 'color', clrs(3,:), 'linewidth', 2.5)

if ifAddGaussianFits
 pars0 = [20;0];
 fitoptions = [];
 fitoptions.DerivativeCheck = 'on';
 fun = @(pars) computeGradient2oMaxEnt(pars, n, mu1, mu2);
 [parsOut,~] = minFunc(fun, pars0, fitoptions);
 [f,g,Ex,Ex2] = computeGradient2oMaxEnt(parsOut,n, mu1, mu2);
 disp('desired and true E(K) and E(K^2)')
  [mu1, Ex; mu2, Ex2],  pause; 
 pcount_gauss = exp(parsOut(1) * (0:n) + parsOut(2) * (0:n).^2);
 pcount_gauss = pcount_gauss/sum(pcount_gauss);
 semilogy(0:n, pcount_gauss, ':', 'color', clrs(3,:), 'linewidth', 2.5)
end
clear data_1 data_2 mu1 mu2 a b Z lognchoosek logpcount ks pcount

% last, add data from Tkacik et al., 2012 
% ('The simplest maximum entropy model for...')
n = 40;
data_1=load('../data/other_studies/figure_data/Tkacik_Real_Trace.mat');
data_1=data_1.Tkacik_real_trace;
ks = round(data_1(:,1));
pcount = data_1(:,2); pcount = pcount/sum(pcount);
mu1 =  (ks') * pcount;               mu2 = (ks').^2 * pcount;
Z = ( n * (mu2/mu1 - mu1 -1)) + mu1;
a = (n * mu1 - mu2) / Z;             b = (n - mu1) * (n - mu2/mu1) / Z;
lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';
logpcount = lognchoosek + betaln(a + (0:n), n + b - (0:n))' - betaln(a, b);
[a2,b2] = fitBetabinML(pcount, n, 100000);
logpcount2 = lognchoosek + betaln(a2 + (0:n), n + b2 - (0:n))' ...
             - betaln(a2, b2);

semilogy(data_1(:,1), pcount, 's-', 'color', clrs(5,:), ...
         'linewidth', 1.5, 'markerSize', 1.5)
hold on
semilogy(data_1(:,1), pcount, 's', 'color', clrs(5,:), ...
         'linewidth', 0.75, 'markerSize', 3)
semilogy(0:n, exp(logpcount), '-', 'color', clrs(5,:), 'linewidth', 2.5)

if ifAddGaussianFits
 pars0 = [20;0];
 fitoptions = [];
 fitoptions.DerivativeCheck = 'on';
 fun = @(pars) computeGradient2oMaxEnt(pars, n, mu1, mu2);
 [parsOut,~] = minFunc(fun, pars0, fitoptions);
 [f,g,Ex,Ex2] = computeGradient2oMaxEnt(parsOut,n, mu1, mu2);
 disp('desired and true E(K) and E(K^2)')
  [mu1, Ex; mu2, Ex2],  pause; 
 pcount_gauss = exp(parsOut(1) * (0:n) + parsOut(2) * (0:n).^2);
 pcount_gauss = pcount_gauss/sum(pcount_gauss);
semilogy(0:n, pcount_gauss, ':', 'color', clrs(5,:), 'linewidth', 2.5)
end
clear data_1 data_2 mu1 mu2 a b Z lognchoosek logpcount ks pcount

set(gca, 'FontSize', fontSize)  
xlabel('population spike count', 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight )
axis([0, 61, 1.1*10^(-4), 1]); 
box off, set(gca, 'TickDir' ,'out')

%% c) goodness of beta-binomial fit
load('fig_data/fig2_data_alphabeta.mat')

if splitFigs
  figure(24)
  subplot(1,2,1)
else
  subplot(21,34,vec(bsxfun(@plus,  (11:18)', (0:4)*34)))
end
h = area(Ns, [mean(a,2)' - sqrt(var(a')); 2*sqrt(var(a'))]'); 
h(2).FaceColor = clrs(1,:); h(1).FaceColor = [1,1,1];
h(2).EdgeColor = [1,1,1]; h(1).EdgeColor = [1,1,1];
hold on
plot(Ns, mean(a,2), 'color', clrs(5,:), 'linewidth', 3),
hold off
set(gca, 'TickDir', 'out')
set(gca, 'Xtick', [])
set(gca, 'Ytick', [0.35,0.4,0.45])
box off
set(gca, 'Linewidth', axesThickness)
axis([Ns(1), Ns(end), 0.32,0.49])
set(gca, 'FontSize', fontSize)
ylabel('\alpha', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )

if splitFigs
  figure(24)
  subplot(1,2,2)
else
  subplot(21,34,vec(bsxfun(@plus,  (11:18)', (5:9)*34)))
end
h = area(Ns, [mean(b,2)' - sqrt(var(b')); 2*sqrt(var(b'))]'); 
h(2).FaceColor = clrs(1,:); h(1).FaceColor = [1,1,1];
h(2).EdgeColor = [1,1,1]; h(1).EdgeColor = [1,1,1];
hold on
plot(Ns, mean(b,2), 'color', clrs(5,:), 'linewidth', 3)
set(gca, 'TickDir', 'out')
set(gca, 'Xtick', [20, 160, 300])
set(gca, 'Ytick', [11,13,15])
box off
set(gca, 'Linewidth', axesThickness)
axis([Ns(1), Ns(end), 10.8,15.2])
set(gca, 'FontSize', fontSize)
xlabel('population size','FontName',fontName,'FontSize',fontSizeXlabel, ...
    'FontWeight', fontWeight )
ylabel('\beta', 'FontName',fontName,'FontSize',fontSizeYlabel,...
    'FontWeight', fontWeight )


%% d) Heat traces flat model


if splitFigs
  figure(23)
else
  subplot(21,34,vec(bsxfun(@plus,  (21:34)', (0:15)*34)))
end

load('../results/K_pairwise_final/heat_traces_flat_nat.mat')
maxi = max(ns/10)/2;
c = cN(2:2:end,:,:);

for i = 1:maxi,
    plot(Ts, squeeze(mean(c(i,:,:),2)),'color',clrs(i,:),'linewidth',2.5);
    hold on;
end
for i = 1:maxi,
    plot(Ts, squeeze(c(i,:,:)), '-', 'color', clrs(i,:), 'linewidth', 0.8);
    hold on;
end

line([1,1], [0, 1.05*max(c(:))], 'linestyle', ':', 'color', 'k', ...
    'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([min(Ts),max(Ts),0.95*min(c(:)), 1.05*max(c(:))]); %axis autoy
set(gca, 'XTick', [1, 1.5, 2]);
set(gca, 'YTick', [1,3,5])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('temperature','FontName',fontName,'FontSize',fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('specific heat c','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )

%% inset
inset2_1  = figure('Tag','fig2+1','units','centimeters','position',...
    [0,0,4,3]);
maxi = size(cN,1)/2; [~,idxT] = min(abs(Ts-1));
lclrs = [0*[1,1,1]; 0.5*[1,1,1]];
     
cms = max(c(:,:,:),[],3);
c1s = c(:,:,idxT) + (1-Ts(idxT))/(Ts(idxT+1)-Ts(idxT)) * ...
    (c(:,:,idxT+1)-c(:,:,idxT));
plot(20, mean(cms(1,:)), 'color', lclrs(1,:), 'linewidth', 1.5);
hold on
plot(20, mean(c1s(1,:)), 'color', lclrs(2,:), 'linewidth', 1.5);
     
% plot max{c(T)}
plot(20*(1:maxi), mean(cms,2), '-', 'color', lclrs(1,:), 'linewidth', 1.5)
for i = 1:maxi, plot(20*i, mean(cms(i,:),2), 's','color', clrs(i,:), ...
        'linewidth', 2, 'markerSize', 2); end

% plot c(T=1)
plot(20*(1:maxi), mean(c1s,2), '-', 'color',lclrs(2,:), 'linewidth', 1.5)
for i = 1:maxi, plot(20*i, mean(c1s(i,:),2), 's','color', clrs(i,:), ...
        'linewidth', 2, 'markerSize', 2); end

set(gca, 'TickDir', 'out')
set(gca, 'Xtick', [20, 60, 100])
set(gca, 'Ytick', [2, 4])
box off
set(gca, 'Linewidth', axesThickness)
axis([15,20*maxi+5, 0.45, 5])
set(gca, 'FontSize', fontSize)
xlabel('n', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
ylabel('c', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )
legend('max c', 'c at T=1', 'location', 'Northwest') 
legend boxoff

hold off


%% e) Heat traces binomial and flat ising


if splitFigs
  figure(25)
  subplot(1,2,1)
else
  figure(figure2)
  subplot(21,34,vec(bsxfun(@plus,  (21:26)', (18:20)*34)))
end

load('fig_data/fig2_data.mat')
maxi = length(Ns);

for i = 1:maxi,
  plot(Ts, squeeze(mean(cN(i,:,:,1),3)),'color',clrs(i,:),'linewidth',2.5);
  hold on;
end
for i = 1:maxi,
    plot(Ts, squeeze(cN(i,:,:,1)), '-', 'color',clrs(i,:),'linewidth',0.8);
    hold on;
end
line([1,1], [0, 1.05*max(cN(:))], 'linestyle', ':', 'color', 'k', ...
'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([min(Ts),max(Ts),0.95*min(vec(squeeze(cN(:,:,:,1)))), ...
    1.05*max(vec(squeeze(cN(:,:,:,1))))]); %axis autoy
set(gca, 'XTick', [1, 1.5, 2]);
set(gca, 'YTick', [0.3,0.4])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('temperature','FontName',fontName,'FontSize',fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('specific heat','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )


if splitFigs
  figure(25)
  subplot(1,2,2)
else
  subplot(21,34,vec(bsxfun(@plus,  (29:34)', (18:20)*34)))
end

for i = 1:maxi,
  plot(Ts, squeeze(mean(cN(i,:,:,2),3)),'color',clrs(i,:),'linewidth',2.5);
  hold on;
end
for i = 1:maxi,
    plot(Ts, squeeze(cN(i,:,:,2)),'-','color',clrs(i,:),'linewidth', 0.8);
    hold on;
end

line([1,1], [0, 1.05*max(cN(:))], 'linestyle', ':', 'color', 'k', ...
    'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([min(Ts),max(Ts),0.95*min(vec(squeeze(cN(:,:,:,2)))), ...
    1.05*max(vec(squeeze(cN(:,:,:,2))))]); %axis autoy
set(gca, 'XTick', [1, 1.5, 2]);
set(gca, 'YTick', [0.4,0.8,1.2])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('temperature','FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('specific heat', 'FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )


%% suppl. figure recently moved to main paper, figure 3 (d-f):

Tmin = 0.5; 
Tmax = 3.0;

Ns = [100];                           % population sizes
qs = [0.5,1,1.5,2,2.5,3,3.5,4, ...
      50*0.0832, ...
      4.5,5,5.5,6,6.5,7,7.5,8,8.5,9]/50; % firing probabilities
rhos = [0,0.00001,0.000025,0.00005,0.000075,0.0001,0.00025,0.0005,...
        0.001,0.00175,0.0025,0.00375,0.005,0.00675,0.0075,0.00825,...
        0.01,0.015,0.02,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25];
Tinv = 1./linspace(Tmin, Tmax, 1000);   % inverse temperatures


res = zeros(numel(qs),numel(rhos));
loc = zeros(numel(qs),numel(rhos));

ifPlot = 1;

transitions = cell(length(Ns),1);
locations   = cell(length(Ns),1);
options = optimoptions('fminunc');
options.Algorithm = 'quasi-newton';
options.TolFun = 1e-10;
options.TolX = 1e-10;
options.Display = 'off';
if ifPlot
    figure;
end    
   
for ni = 1:length(Ns)
    n = Ns(ni);
    lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))'; 
    res = zeros(numel(qs),numel(rhos));
    
    for j = length(rhos):-1:1
      disp(rhos(j))
      for i = 1:length(qs)
        rho = rhos(j);
        q = qs(i);

        if rho==0
            % independent case: statistics are independent of n, i.e. we can
            % look at a Binomial distribution with n=1 (-> Bernoulli), varying T 
            qb = q.^Tinv ./ (q.^Tinv + (1-q).^Tinv);    
            c = qb.* (1-qb) .* (log(qb) - log(1-qb)).^2;
            error = 2 * ((1-q).^Tinv + q.^Tinv) - log(q/(1-q)) * Tinv .* (q.^Tinv - (1-q).^Tinv); 
            %disp(' rho = 0 ')
        else
            % correlated case: use beta-binomial for fixed n=100. 
            %
            mu1 =  n * q; % (0:n) * pcount;
            mu2 =  n * q * (n * q + (1-q) * (1 + (n-1) * rho)); %(0:n).^2 * pcount;
            Z = ( n * (mu2/mu1 - mu1 -1)) + mu1;
            a = (n * mu1 - mu2) / Z;
            b = (n - mu1) * (n - mu2/mu1) / Z;    
            logpcount = lognchoosek + betaln(a + (0:n), n+b-(0:n))'-betaln(a,b);
            pcount = exp(logpcount);    
            [c, pcountT] = computeSpecificHeatTraces(pcount, 1./Tinv, 0);
            %disp( mu2 - mu1^2 )

        end

        %disp( (0:n).^2 * pcount - ((0:n) * pcount).^2 )
        
        f = @(z) - interp1(1./Tinv, c, z, 'spline');
        [~, optb] = max(c); 
        res(i,j) = fminunc(f, 1./Tinv(optb), options);
        loc(i,j) = 1/Tinv(optb);
        if 1./Tinv(optb) == res(i,j) 
            disp('warning: fminunc didnt do anything')
        end
        
        if ifPlot
            subplot(3,1,2)
            plot(1./Tinv, c)
            hold on
            plot(1./Tinv(optb), c(optb), 'g*')        
            line([1, 1], [0, 1.1*max(c)], 'linewidth', 2, 'color', 'k')
            title('c(T)')
            xlabel(['firing rate = ', num2str(q*50), ' Hz'])
            axis([Tmin, Tmax, 0.2, 1.1*max(c)])
            hold off
            if rho == 0
                subplot(3,1,3)
                plot(1./Tinv, abs(error))
                xlabel('T')
                title('minimum of error')
                axis([Tmin, Tmax, 0, 1.1*max(abs(error))])
                subplot(3,1,2)
                hold on
                plot(1./Tinv(optb), qb(optb).* (1-qb(optb)) .* (log(qb(optb)) - log(1-qb(optb))).^2, 'r*')
                hold off
            end
            pause(0.1);        
        end 
      end
    end
    transitions{ni} = res;
    locations{ni} = loc;
end   
qs = 50*qs;

figure;
subplot(1,2,1)
log0 = 10^(-6); % log(0) for plotting purposes

qs_idx_show = [4, 9, 13];
clrs = [254,153,41;
236,112,20;
204, 76,2;      % colors for inidivual traces.
153,52,4;
102,37,6;
0,0,0]/255;


for i = 1:length(qs_idx_show)
    plot(0.0001, 0.0001, 'linewidth', 2, 'color', clrs(i,:))
    hold on
end
line([log0, 0.25], [1,1], 'color', 'k', 'linewidth', 1.5)
lgnd = cell(size(qs_idx_show));
for i = 1:length(qs_idx_show), 
    plot([rhos(1)+log0,rhos(2)], ...
             [locations{1}(qs_idx_show(i),1), locations{1}(qs_idx_show(i),2)], ...
             '--', 'linewidth', 2, 'color', clrs(i,:)); 
    plot(rhos(2:end), locations{1}(qs_idx_show(i),2:end), ...
        'linewidth', 2, 'color', clrs(i,:)); 
    lgnd{i} = [num2str(qs(qs_idx_show( i ))), 'Hz'];
end
legend(lgnd)
xlabel('\rho'); 
ylabel('T*'); 
title('T*(\rho) for fixed FR'); 
axis([-0.005, 0.25, 0.8, 1.35])
box off
set(gca, 'TickDir', 'out')
legend('boxoff')
%set(gca, 'XTick', [10^(-5), 10^(-3), 10^(-1)])

subplot(1,2,2)
rho_idx_show = [1, 17, 20, 21, 23, 28];
clrs = flipud(copper(2*length(rho_idx_show)));
for i = 1:length(rho_idx_show)
    plot(0.0001, 0.0001, 'linewidth', 2, 'color', clrs(i,:))
    hold on
end
line([0, 9], [1,1], 'color', 'k', 'linewidth', 2)
lgnd = cell(size(rho_idx_show));
for i = 1:length(rho_idx_show),
    
    plot(qs, locations{1}(:,rho_idx_show(i)),  ...
         'linewidth', 2, 'color', clrs(i,:)); 
    hold on; 
    lgnd{i} = [num2str(rhos(rho_idx_show( i )))];
end; 
legend(lgnd)

title('T*(FR) for fixed \rho'); 
xlabel('FR'); 
ylabel('T*'); 
line([4.16,4.16], [0.6, 1.65], 'color', 'k', 'lineStyle', '--')
axis([0, 9, 0.8, 1.35]);
axis autoy
box off
set(gca, 'TickDir', 'out')
legend('boxoff')




load('fig_data/fig3_peakHeatData.mat')
clrs = [254,153,41;236,112,20;204,76,2;153,52,4;102,37,6;0,0,0]/255;
options = optimoptions('fminunc');
options.Algorithm = 'quasi-newton';
options.TolFun = 1e-10;
options.TolX = 1e-10;
options.Display = 'off';

figure
for ni = 1:length(Ns)
    res = transitions{ni};
    n = Ns(ni);
    qopts = zeros(length(rhos),1);
    for j = 1:length(rhos)
        f =  @(z) abs( interp1(qs(1:end-1), res(1:end-1,j), z) - 1);
        [~, optq] = min( abs( res(1:end-1,j) - 1)); 
        qopts(j) = fminunc(f, qs(optq), options);
    end
    res_plot = 2*ones(numel(qs)-1,numel(rhos));
    idx = res(1:end-1,:)==1;
    res_plot(idx) = 1;
    idx = res(1:end-1,:)>1;
    res_plot(idx) = 0;       
    
    plot(qopts, rhos, 'w', 'linewidth', 2, 'color', clrs(ni,:))
    hold on
    
end
xlabel('firing rate [Hz]')
ylabel('correlation \rho')
line([4.16, 4.16], [0, 10], 'linestyle', '--', 'color', 'k')
axis([0, 9, 0, 0.25])
set(gca, 'XTick', 0:4:8)
set(gca, 'YTick', 0:0.1:0.2)
box off
set(gca, 'TickDir', 'out')
legend('n = 20', 'n = 40', 'n = 60', 'n = 80', 'n = 100', 'n = 120', 'location', 'NorthWest') 


%% export figure to specified size
% if ~splitFigs
%              s.Version= '1';
%              s.Format= 'pdf';
%             s.Preview= 'none';
%               s.Width= '19';
%              s.Height= '14';
%               s.Units= 'centimeters';
%               s.Color= 'rgb';
%          s.Background= 'w';
%  %     s.FixedFontSize= 'auto'
%  %    s.ScaledFontSize= 'auto'
%  %          s.FontMode= 'scaled'
%  %       s.FontSizeMin= '8'
%  %    s.FixedLineWidth= '1'
%  %   s.ScaledLineWidth= 'auto'
%  %          s.LineMode= 'scaled'
%  %      s.LineWidthMin= '2'
%  %          s.FontName= 'auto'
%  %        s.FontWeight= 'auto'
%  %         FontAngle= 'auto'
%  %      s.FontEncoding= 'auto'
%  %           s.PSLevel= '2'
%  %          s.Renderer= 'auto'
%  %        s.Resolution= 'auto'
%  %      s.LineStyleMap= 'none'
%          s.ApplyStyle= '1';
%              s.Bounds= 'loose';
%  %          s.LockAxes= 'on'
%  %            s.ShowUI= 'on'
%  %      s.SeparateText= 'off'
%        
%  figure(figure2);
%  hgexport(figure2,'f2.pdf',s);     
%  figure(figure2);
%  pause;
%  export_fig('fig2.pdf')
%  
%  
%               s.Width= '4';
%              s.Height= '3';
% 
%  figure(inset2_1);             
%  hgexport(inset2_1,'test.pdf',s);     
%  figure(inset2_1);            
%  pause; 
%  export_fig('fig2_1.pdf')
% % hgexport(inset2_2,'test.pdf',s);     
% % export_fig('fig2_2.pdf')
% % hgexport(inset2_3,'test.pdf',s);     
% % export_fig('fig2_3.pdf')
%  
%  
% end
