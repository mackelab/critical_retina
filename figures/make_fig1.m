%% produces fig 1 for the journal version of 'criticality and correlations' 
% figure '1' is the overview figure summarising our simulation and how
% we find signatures of 'criticality' in this data

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

%%

% subplot b) Plot RGC centers
if splitFigs
  figure(12)
else
  subplot(21,34,vec(bsxfun(@plus,  (1:10)', (14:20)*34)))
%  subplot(17,11,vec(bsxfun(@plus, (7:11)', (0:4)*11)))
end
magFactor = 1.05;
clrs = copper(2*6000);  % choose colors in analogy to 
clrs = clrs(end/2+1:end,:); % figure 4 - RGC center plot
clrs = clrs(3/4*end,:); %
plot(double(pixel_size)*cell_positions(1,:), ...
     double(pixel_size)*cell_positions(2,:), ...
     'o', 'color', clrs, 'linewidth', 1.25, 'markerSize', 3)
axis(magFactor*[-150,150,-100,100])
set(gca, 'YTick', [-50,0,50]);
set(gca, 'XTick', [-100,0,100]);
set(gca, 'XAxisLocation', 'top')
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


%% subplot c) Make raster plot
if splitFigs
  figure(13)
else
  subplot(21,34,vec(bsxfun(@plus, 2*(7:11)', (0:8)*34)))
end
raster = output.spikes(101:200,1000:1100);
imagesc(1-raster);
colormap('gray')
set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick', 1:50:size(raster,2));
set(gca, 'XTickLabel', (0:50:size(raster,2))/50);
YTick = size(raster,1):-50:0;
%YTick = [];
set(gca, 'YTick', YTick(end:-1:1));
set(gca, 'YTickLabel', floor(size(raster,1)/50)*50:-50:0);
set(gca, 'FontSize', fontSize)
axis([-1, size(raster,2), 0.5, size(raster,1)+1.5])
xlabel('time (s)', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
ylabel('neuron label','FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )
box off, set(gca, 'TickDir' ,'out')

%% subplot d) 
% Plot firing rate 
if splitFigs
  figure(14)
  subplot(2,2,1)
else
  subplot(21,34,vec(bsxfun(@plus, (14:16)', (11:14)*34)))
end
binSize = 50; % 20ms bins => factor 1/0.02 = 50 
firingRates = binSize*mean(output.spikes,2);
xrange = linspace(0.95*min(firingRates),1.05*max(firingRates),25);
h = histc(firingRates, xrange);
h = h/sum(h);

clrs = [254,153,41;236,112,20;204, ...     % colors for inidivual traces.
        76,2;153,52,4;102,37,6;0,0,0]/255; % Currently picked by hand.

bar(xrange, h, 1, 'faceColor', clrs(2,:),'edgeColor', clrs(2,:));
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickLength', 2.5*get(gca, 'TickLength'))
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

% Plot correlation coefficient histograms
if splitFigs
  figure(14)
  subplot(2,2,2)
else
  subplot(21,34,vec(bsxfun(@plus, (14:16)', (17:20)*34)))
end
idxM = logical(triu(ones(size(output.spikes,1)),1));
corrs = output.spkCorrs(idxM);
xrange = linspace(min(corrs),1.05*max(corrs),40);
h = histc(corrs, xrange);
h = h/sum(h);
bar(xrange, h, 1, 'faceColor', clrs(2,:), 'edgeColor', clrs(2,:));
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickLength', 2.5*get(gca, 'TickLength'))
if min(corrs) < 0
 mx = 2*min(corrs); 
else
 mx = 0.9*min(corrs); 
end    
Mx = 1*max(corrs);
axis([mx,Mx,-0.05*max(h),1.05*max(h)])
set(gca, 'YTick', [0:0.1:max(h)]);
set(gca, 'XTick', 0:0.2:0.4);
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
txt = ['{\fontname{',fontName,'}\fontsize{', num2str(fontSizeXlabel), ...
    '}\rho_{ij}}'];
xlabel(['correlation'], 'FontName',fontName,'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )

% Plot P(K)
if splitFigs
  figure(14)
  subplot(2,2,3:4)
else
  subplot(21,34,vec(bsxfun(@plus, (18:22)', (11:20)*34)))
end
K = sum(output.spikes,1); 
xrange = 0:max(K)+1; 
h = histc(K, xrange);
h = h/sum(h);
semilogy(xrange, h, '-', 'color', clrs(2,:),'linewidth', 1.5, ...
    'markerSize', 1.5);
hold on
semilogy(xrange, h, 's', 'color', clrs(2,:),'linewidth', 1.5, ...
    'markerSize', 1.5);

set(gca, 'Linewidth', axesThickness)
mx = -1; Mx = 120; %1.05*max(K);
axis([mx,Mx,10^(-5.2),10^(0)]); %axis autoy
set(gca, 'XTick', 0:40:max(K));
set(gca, 'YTick', 10.^[-5,-3,-1])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('population spike count'   , 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('frequency', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )

%% e) googness of fit

% used temperatures (slight annoyance with getting the ordering right)
Ts = [ 0.8000, 0.8563, 0.9000, 0.9363, 0.9675, 0.9938, 1.0175, 1.0388, ...
       1.0388, 1.0575, 1.0775, 1.0988, 1.1200, 1.1413, 1.1650, 1.1900, ...
       1.2150, 1.2425, 1.2713, 1.3013, 1.3338, 1.3687, 1.4063, 1.4463, ...
       1.4463, 1.4900, 1.5375, 1.5900, 1.6500, 1.7175, 1.7950, 1.8875, ...
       2.0000];

n = 100;
idxRepets = [1];
range = (100*Ts(1):100*Ts(end));
clrs = jet( length(range) );

c = zeros(31,1);

if splitFigs
  figure(22)
  subplot(1,2,1)
else
  subplot(21,34,vec(bsxfun(@plus,  (25:31)', (0:4)*34)))
end

for idxRepet = idxRepets 
 for i = 1:31,
  tmp = num2str(Ts(i)); tmp(tmp=='.')='_';    
  load(['../results/method_validation/specific_heat_samples_long_runs/',...
        'longRun',num2str(n),'idxRepet',num2str(idxRepet),'T',tmp,'.mat']) 
  x = linspace( 0.001, (n/100)^2 * 12, size(MoE,2)-1);
  tmp = MoE(2,1:end-1) - MoE(1,1:end-1).^2;
  tmp = cumsum(tmp) ./ (1:length(tmp));
  [~, idxClr] = min(abs( 100*Ts(i) - range));
  semilogx(x, tmp/n, 'color', clrs(idxClr,:), 'linewidth', 0.75),
  hold on
  c(i) = mean(MoE(2,2:end-1)-MoE(1,2:end-1).^2)/n;
 end
end
set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick',      [ 1/64,   1/32,   1/16,   1/8,   ...
                          1/4,   1/2,   1,  2,  4,  8])
set(gca, 'XTickLabel', {'1/64', '1/32', '1/16', '1/8', ...
                         '1/4', '1/2', '1','2','4','8'})
set(gca, 'YTick', [0, 0.5, 1])
set(gca, 'FontSize', fontSize)  
xlabel('time [h]', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
ylabel('specific heat', 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
box off, set(gca, 'TickDir' ,'out')
axis([1/60,12,0,1.45])
hold off

if splitFigs
  figure(22)
  subplot(1,2,2)
else
  subplot(21,34,vec(bsxfun(@plus,  (32:34)', (0:4)*34)))
end
plot(Ts(1:31), c, 'color', 0.3*[1,1,1], 'linewidth', 1.25);
hold on
for idxRepet = idxRepets 
 for i = 1:31,
  tmp = num2str(Ts(i)); tmp(tmp=='.')='_';    
  %load(['../results/method_validation/specific_heat_samples/VarEsn',...
  %num2str(n), 'idxRepet', num2str(idxRepet), 'T', tmp, '.mat']) 
  [~, idxClr] = min(abs( 100*Ts(i) - range));
  plot(Ts(i), c(i,idxRepet), 's', 'color', clrs(idxClr,:), ...
      'linewidth', 1.5, 'markerSize', 1.5)
 end
end
hold off
set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick', [1, 1.5, 2])
set(gca, 'FontSize', fontSize)
xlabel('temperature', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
set(gca, 'YTick', [])
box off, set(gca, 'TickDir' ,'out')
axis([0.8, 2,0,1.45])

%% f) Heat traces
if splitFigs
  figure(16)
else
  subplot(21,34,vec(bsxfun(@plus, (25:34)', (7:20)*34)))
end
maxN = 120;
maxi = maxN / 20;
clrs = [254,153,41;
236,112,20;
204, 76,2;      % colors for inidivual traces.
153,52,4;
102,37,6;
0,0,0]/255;

Ts = [0.8000,0.8563,0.9000,0.9363,0.9675,0.9938,1.0175,1.0388,1.0575, ...
    1.0775,1.0988,1.1200,1.1413,1.1650,1.1900,1.2150,1.2425,1.2713, ...
    1.3013,1.3338,1.3687,1.4063,1.4463,1.4900,1.5375,1.5900,1.6500, ...
    1.7175,1.7950,1.8875,2.0000];

c = cN(2:2:end,:,:);

for i = 1:maxi,
    plot(Ts, squeeze(mean(c(i,:,:),2)),'color',clrs(i,:),'linewidth',2.5);
    hold on;
end
for i = 1:maxi,
    plot(Ts, squeeze(c(i,:,:)), '-', 'color', clrs(i,:), 'linewidth', 0.8);
    hold on;
end
line([1,1], [0, 1.05*max(c(:))], 'linestyle', '--', 'color', 'k', ...
    'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([min(Ts),max(Ts),0.95*min(c(:)), 1.05*max(c(:))]); %axis autoy
set(gca, 'XTick', [1, 1.5, 2]);
set(gca, 'YTick', [0.5, 1, 1.5])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('temperature'   , 'FontName',fontName,'FontSize',fontSizeXlabel, ...
    'FontWeight', fontWeight ) 
ylabel('specific heat c','FontName',fontName,'FontSize',fontSizeYlabel, ...
    'FontWeight', fontWeight )

%% inset
inset1_1=figure('Tag', 'fig1', 'units','centimeters','position',[0,0,4,3]);
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
set(gca, 'Ytick', [0.5, 1,1.5])
box off
set(gca, 'Linewidth', axesThickness)
axis([15,20*maxi+5, 0.45, 1.6])
set(gca, 'FontSize', fontSize)
xlabel('n', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
ylabel('c', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
    'FontWeight', fontWeight )
legend('max c', 'c at T=1', 'location', 'Northwest') 
legend boxoff

hold off

%% export figure to specified size
% if ~splitFigs
%              s.Version= '1';
%              s.Format= 'pdf';
%             s.Preview= 'none';
%               s.Width= '19';
%              s.Height= '11';
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
%  figure(figure1);
%  hgexport(figure1,'/home/mackelab/Desktop/test.pdf',s);     
%  figure(figure1);
%  pause;
%  export_fig('/home/mackelab/Desktop/Projects/fig1.pdf')
%  
% %%
%               s.Width= '4';
%              s.Height= '3';
%              
%  figure(inset1_1);             
%  hgexport(inset1_1,'/home/mackelab/Desktop/test2.pdf',s);     
%  figure(inset1_1);             
%  pause;
%  export_fig('/home/mackelab/Desktop/Projects/fig1_1.pdf')
% 
%  
% end
