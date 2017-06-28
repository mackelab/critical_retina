%% produces fig 1 for the journal version of 'criticality and correlations' 
% figure '1' is the overview figure summarising the criticality analysis
% procedure introduced by Tkacik et al. 

clear all

splitFigs = false;

if ~splitFigs
  figure1 = figure('Tag', 'fig1', 'units','centimeters','position',...
      [0,0,19,11]);
end

axesThickness  = 1.1; % overall thickness of axes. Doesn't seem to do much.

fontName = 'Arial';    fontWeight     = 'normal';
fontSize       = 1 * 10;   fontSizeTitle  = 1 * 16;   
fontSizeXlabel = 1 * 10;   fontSizeYlabel = 1 * 11;
fontSizeText   = 1 * 10;   fontSizeLegend = 1 * 11;


% Load data and preprocess. What is needed after this is a struct with 
% spike raster output.spikes (n x t matrix, with t being total recording
% time), and spike correlations and covariance in output.spkCorrs and
% output.spkCovs
load('../data/RGC_sim_nat.mat')        % loads raw simulated RGC data
load('../results/K_pairwise_final/heat_traces_nat.mat') % specific heat 
load('../results/K_pairwise_final/idxSubsamples.mat')

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

% pick two population sizes: n1, n2. 
% As figure 1 serves mostly just as an illustration, also pick the example
% subpopulations per hand: The selection below was made because the 
% histograms over firing rates have similar peak values and can thus be 
% plotted with the same bar hights and y-axis. 
n1 = 2; 
n2 = 6;
i1 = 7;
i2 = 1;

%%

% subplot a) top: Plot RGC centers
if splitFigs
  figure(12)
else
   subplot(2,4,1)

end
magFactor = 1.05;

clrs = [254,153,41;236,112,20;204, ...     % colors for inidivual traces.
        76,2;153,52,4;102,37,6;0,0,0]/255; % 

plot(double(pixel_size)*cell_positions(1,:), ...
     double(pixel_size)*cell_positions(2,:), ...
     'o', 'color', 0.5*[1,1,1], 'markerSize', 3, ...
     'markerFaceColor', 0.5*[1,1,1])
hold on 

idx1 = idxSubsamples{n1}(:,i1);
plot(double(pixel_size)*cell_positions(1,idx1), ...
     double(pixel_size)*cell_positions(2,idx1), ...
     'o', 'color', clrs(n1/2,:), 'markerSize', 3, ...
     'markerFaceColor', clrs(n1/2,:))
idx2 = idxSubsamples{n2}(:,i2);
plot(double(pixel_size)*cell_positions(1,idx2), ...
     double(pixel_size)*cell_positions(2,idx2), ...
     'o', 'color', clrs(n2/2,:), 'markerSize', 3, ...
     'markerFaceColor', clrs(n2/2,:))
hold off
axis(magFactor*[-150,150,-100,100])
set(gca, 'YTick', [-50,0,50]);
set(gca, 'XTick', [-100,0,100]);
set(gca, 'XAxisLocation', 'bottom')
set(gca, 'visible', 'on')
set(gca, 'Linewidth', axesThickness)
box off; set(gca, 'TickDir' ,'out')
line(magFactor*[-150,150], magFactor*[100,100], 'color', 'k', ...
    'linewidth', axesThickness)
line(magFactor*[-150,150], -magFactor*[100,100], 'color', 'k', ...
    'linewidth', axesThickness)
line(magFactor*[150,150], magFactor*[-100,100], 'color', 'k', ...
    'linewidth', axesThickness)
line(-magFactor*[150,150], magFactor*[-100,100], 'color', 'k', ...
    'linewidth', axesThickness)
txt=['{\fontname{',fontName,'}\fontsize{',num2str(fontSizeXlabel),'}\mu}'];
set(gca, 'FontSize', fontSize)
xlabel(['x-coordinate [', txt,'m]'], 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
ylabel(['y-coordinate [', txt,'m]'], 'FontName', fontName, ...
    'FontSize', fontSizeYlabel, 'FontWeight', fontWeight )

title('subsample population recordings','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )


%% subplot a, bottom) Make raster plot
if splitFigs
  figure(13)
else
  subplot(2,4,5)
  %subplot(21,34,vec(bsxfun(@plus, 2*(7:11)', (0:8)*34)))
end
idx1 = sort(idxSubsamples{n1}(:,i1));
idx2 = sort(idxSubsamples{n2}(:,i2));

idx12 = sort(randsample(n, 100, false));
raster = output.spikes(idx12,1000:1050);

for i = 1:size(raster,1)
    if ismember(i, idx1)
        clr = clrs(n1/2,:); 
    elseif ismember(i, idx2),
        clr = clrs(n2/2,:); 
    else
        clr = 0.5 * [1,1,1];
    end
    plot(find(raster(i,:)), i * ones(sum(raster(i,:)>0),1), 's', ...
    'color', clr, 'markerFaceColor', clr, 'markersize', 4, 'linewidth', 0.5);
    hold on
end
colormap([clrs(n1/2,:); clrs(n2/2,:);[1,1,1]])
axis([-1, size(raster,2), -0.5, size(raster,1)+1.5])

set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick', 1:50:size(raster,2));
set(gca, 'XTickLabel', (0:50:size(raster,2))/50);
YTick = [100;1];
set(gca, 'YTick', size(raster,1)-YTick);
set(gca, 'YTickLabel', {'100', '1'});
set(gca, 'FontSize', fontSize)
xlabel('time (s)', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
ylabel('neuron label','FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )
box off, set(gca, 'TickDir' ,'out')


%% subplot d) 
% Plot firing rate 
binSize = 50; % 20ms bins => factor 1/0.02 = 50 
firingRates = binSize*mean(output.spikes,2);
xrange = linspace(0.95*min(firingRates),1.05*max(firingRates),25);
h = histc(firingRates, xrange);
sh = sum(h);
h = h/sum(h);

idx1 = idxSubsamples{n1}(:,i1);
h1 = histc(firingRates(idx1), xrange);
sh = sum(h1);
h1 = h1/sh;
idx2 = idxSubsamples{n2}(:,i2);
h2 = histc(firingRates(idx2), xrange);
sh = sum(h2);
h2 = h2/sh;

clrs = [254,153,41;236,112,20;204, ...     % colors for inidivual traces.
        76,2;153,52,4;102,37,6;0,0,0]/255; % 
    


subplot(3,8,3)
bar(xrange, h1, 1, 'faceColor', clrs(n1/2,:),'edgeColor', clrs(n1/2,:));
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickLength', 2.5*get(gca, 'TickLength'))
axis square
if min(firingRates) < 0
 mx = 1.2*min(firingRates); 
else
 mx = 0.5*min(firingRates); 
end    
Mx = 1.05*max(firingRates);
axis([mx,Mx,-0.05*max(h),1.05*max(h)])
set(gca, 'YTick', [0:0.1:max(h)]);
set(gca, 'XTick', [1:4]);
set(gca, 'Linewidth', axesThickness)
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('firing rate [Hz]',    'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )
text(3, 0.15, ['n = ', num2str(10*n1)])

title('extract statistics','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )


subplot(3,8,4)
bar(xrange, h2, 1, 'faceColor', clrs(n2/2,:),'edgeColor', clrs(n2/2,:));
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickLength', 2.5*get(gca, 'TickLength'))
axis square
if min(firingRates) < 0
 mx = 1.2*min(firingRates); 
else
 mx = 0.5*min(firingRates); 
end    
Mx = 1.05*max(firingRates);
axis([mx,Mx,-0.05*max(h),1.05*max(h)])
set(gca, 'YTick', [0:0.1:max(h)]);
set(gca, 'XTick', [1:4]);
set(gca, 'Linewidth', axesThickness)
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('firing rate [Hz]',    'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
text(3, 0.15, ['n = ', num2str(10*n2)])


% Plot correlation coefficient histograms
idxM = logical(triu(ones(size(output.spikes,1)),1));
corrs = output.spkCorrs(idxM);
xrange = linspace(min(corrs),1.05*max(corrs),40);
h = histc(corrs, xrange);
sh = sum(h);
h = h/sum(h);

idx1 = zeros( size(output.spikes,1) );
idx1(idxSubsamples{n1}(:,i1)) = 1;
idx1= logical(triu(idx1*idx1',1));
h1 = histc(output.spkCorrs(idx1), xrange);
sh = sum(h1);
h1 = h1/sh;

idx2 = zeros( size(output.spikes,1) );
idx2(idxSubsamples{n2}(:,i2)) = 1;
idx2 = logical(triu(idx2*idx2',1));
h2 = histc(output.spkCorrs(idx2), xrange);
sh = sum(h2);
h2 = h2/sh;

subplot(3,8,11)
bar(xrange, h1, 1, 'faceColor', clrs(n1/2,:),'edgeColor', clrs(n1/2,:));

axis square
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickLength', 2.5*get(gca, 'TickLength'))
if min(corrs) < 0
 mx = 2*min(corrs); 
else
 mx = 0.9*min(corrs); 
end    
Mx = 1*max(corrs);
axis([mx,Mx,-0.05*max(h1),1.05*max(h1)])
set(gca, 'YTick', [0:0.1:max(h1)]);
set(gca, 'XTick', 0:0.2:0.4);
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
txt = ['{\fontname{',fontName,'}\fontsize{', num2str(fontSizeXlabel), ...
    '}\rho_{ij}}'];
xlabel(['correlation'], 'FontName',fontName,'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )

subplot(3,8,12)
bar(xrange, h2, 1, 'faceColor', clrs(n2/2,:),'edgeColor', clrs(n2/2,:));

axis square
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickLength', 2.5*get(gca, 'TickLength'))
if min(corrs) < 0
 mx = 2*min(corrs); 
else
 mx = 0.9*min(corrs); 
end    
Mx = 1*max(corrs);
axis([mx,Mx,-0.05*max(h1),1.05*max(h1)])
set(gca, 'YTick', [0:0.1:max(h2)]);
set(gca, 'XTick', 0:0.2:0.4);
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
txt = ['{\fontname{',fontName,'}\fontsize{', num2str(fontSizeXlabel), ...
    '}\rho_{ij}}'];
xlabel(['correlation'], 'FontName',fontName,'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight ) 


% Plot P(K)
K = sum(output.spikes,1); 
xrange = 0:max(K)+1; 
h = histc(K, xrange);
sh = sum(h);
h = h/sum(h);

idx1 = idxSubsamples{n1}(:,i1);
K1 = sum(output.spikes(idx1,:),1);
h1 = histc(K1, xrange);
sh = sum(h1);
h1 = h1/sh;
idx2 = idxSubsamples{n2}(:,i2);
K2 = sum(output.spikes(idx2,:),1);
h2 = histc(K2, xrange);
sh = sum(h2);
h2 = h2/sh;

subplot(3,8,19)
semilogy(xrange, h1, '-', 'color', clrs(n1/2,:),'linewidth', 1.5, ...
    'markerSize', 1.5);
hold on
semilogy(xrange, h1, 's', 'color', clrs(n1/2,:),'linewidth', 1.5, ...
    'markerSize', 1.5);
hold off
axis square
set(gca, 'Linewidth', axesThickness)
mx = -1; Mx = 120; %1.05*max(K);
axis([mx,11,10^(-4.4),1.01]); 
set(gca, 'XTick', 0:5:max(K1));
set(gca, 'YTick', 10.^[-4,-2,0])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('pop. spike count'   , 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )


subplot(3,8,20)
semilogy(xrange, h2, '-', 'color', clrs(n2/2,:),'linewidth', 1.5, ...
    'markerSize', 1.5);
hold on
semilogy(xrange, h2, 's', 'color', clrs(n2/2,:),'linewidth', 1.5, ...
    'markerSize', 1.5);
hold off
axis square
set(gca, 'Linewidth', axesThickness)
mx = -1; Mx = 120; %1.05*max(K);
axis([mx,33,10^(-4.4),1.01]); 
set(gca, 'XTick', 0:10:max(K2));
set(gca, 'YTick', 10.^[-4,-2,0])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('pop. spike count'   , 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 


%% plot empirical vs. model data moments

load('/home/mackelab/Desktop/Projects/Criticality/code/project_current/results/K_pairwise_final/lambda_nat.mat')
load('/home/mackelab/Desktop/Projects/Criticality/results/EfxSimData.mat')
load('/home/mackelab/Desktop/Projects/Criticality/code/project_current/figures/fig_data/fig0_data.mat')

% make sure we use our color-scheme again
clrs = [254,153,41;236,112,20;204,76,2;153,52,4;102,37,6;0,0,0]/255;
for j = 2:-1:1,
    
    n = 10 * ns(j);
    i = is(j);    
    
    fDescrJ = nchoosek(1:n,2)'; 
    covsx = Efx{n/10}((n+1):(n*(n+1)/2),i) - ... % covariance as computed
           (Efx{n/10}(fDescrJ(1, :),i).* ...
            Efx{n/10}(fDescrJ(2, :),i)); % from Efx    for j = js, 
    varsx = Efx{n/10}(1:n,i) .* (1 - Efx{n/10}(1:n,i));
    corsx = covsx ./(sqrt(varsx(fDescrJ(1,:))).*sqrt(varsx(fDescrJ(2,:))));
    Efy = xSampled{j};
    covsy = Efy((n+1):(n*(n+1)/2)) - ... % covariance as computed
           (Efy(fDescrJ(1, :)).* Efy(fDescrJ(2, :))); % from Efx    
    varsy = Efy(1:n) .* (1 - Efy(1:n));
    corsy = covsy ./(sqrt(varsy(fDescrJ(1,:))).*sqrt(varsy(fDescrJ(2,:))));
    
    subplot(3,8,5)
    if j == 2
        line([0, 5], [0, 5], 'color', 'k')
        hold on
        title('fit maximum entropy models','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )
    end
    plot(50 * Efx{n/10}(1:n,i), 50 * Efy(1:n), 'o', 'color', clrs(n/20,:), ...
        'markersize', 3, 'markerFaceColor', clrs(n/20,:))
    axis square
    hold on;
    axis([0.015, 0.08, 0.015, 0.08] * 50)
    set(gca, 'XTick', [1:4])
    set(gca, 'YTick', [1:4])
    set(gca, 'TickDir', 'out')
    xlabel('firing rate [Hz]',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )
    ylabel('firing rate [Hz]',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )
    box off
    set(gca, 'FontSize', fontSize)
    
    
    subplot(3,8,13)
    if j == 2
        line([-0.1, 0.6], [-0.1, 0.6], 'color', 'k')
        hold on
    end
    plot(corsx, corsy, 'o', 'color', clrs(n/20,:), ..., 
        'markersize', 3, 'markerFaceColor', clrs(n/20,:))
    axis square
    axis([-0.05, 0.6, -0.05, 0.6])
    hold on;
    set(gca, 'XTick', [0, 0.2, 0.4])
    set(gca, 'YTick', [0, 0.2, 0.4])
    set(gca, 'TickDir', 'out')
    xlabel('correlation',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )
    ylabel('correlation',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )
    box off
    set(gca, 'FontSize', fontSize)
    
    
    subplot(3,8,21)
    semilogy(0:n, Efx{n/10}(end-n:end,i), '.-', 'color', 0.7 * [1,1,1])
    hold on
    semilogy(0:n, Efy(end-n:end), '.-', 'color', clrs(n/20,:));
    axis square
    axis([0, 33, 10^-4.4, 1.01])
    hold on;
    xlabel('pop. spike-count',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )
    ylabel('probability',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )
    box off
    set(gca, 'FontSize', fontSize)
    set(gca, 'XTick', [0, 10, 20, 30])
    set(gca, 'YTick', [10^-4, 10^-2, 10^0])
    set(gca, 'TickDir', 'out')
    
end


%% f) Heat traces
if splitFigs
  figure(16)
else
  subplot(1,3,3)
end
maxN = 120;
maxi = maxN / 20;
clrs = [254,153,41;
236,112,20;
204, 76,2;      % colors for inidivual traces.
153,52,4;
102,37,6;
0,0,0]/255;

c = cN(2:2:end,:,:);

plot(Ts, squeeze(c(n1/2,i1,:)), '-', 'color', clrs(n1/2,:), 'linewidth', 3.0);
hold on
plot(Ts, squeeze(c(n2/2,i2,:)), '-', 'color', clrs(n2/2,:), 'linewidth', 3.0);
hold off

line([1,1], [0, 1.05*max(c(:))], 'linestyle', '--', 'color', 'k', ...
    'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([min(Ts),max(Ts),0.95*min(c(:)), 1.2]); %axis autoy
set(gca, 'XTick', [1, 1.5, 2]);
set(gca, 'YTick', [0.5, 1, 1.5])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('temperature'   , 'FontName',fontName,'FontSize',fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('specific heat','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )
title('compute specific heat','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )
