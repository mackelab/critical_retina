%% produces figure three for the journal version of our 'criticality and correlations' project
% figure '3' is the figure summarising dependence of specific heat in flat
% models on the average correlation strength within the data, and how
% we can also see this for the K-pairwise model when varying stimuli.

clear all

splitFigs = false;

if ~splitFigs
  figure3 = figure('Tag', 'fig3', 'units','normalized','position',[0,0,0.99,0.99]);
end

addpath(genpath('../code/'))

clrs = [254,153,41;236,112,20;204,76,2;153,52,4;102,37,6;0,0,0]/255;

axesThickness  = 1.1; % overall thickness of axes. Doesn't seem to do much.

fontName = 'Arial';    fontWeight     = 'normal';
fontSize       = 1 * 10;   fontSizeTitle  = 1 * 16;   
fontSizeXlabel = 1 * 10;   fontSizeYlabel = 1 * 11;
fontSizeText   = 1 * 10;   fontSizeLegend = 1 * 11;

% subplot a) Heat traces for various correlation strengths rho
greyLevel = 0.8;

load('fig_data/fig3_data.mat')
idxSP = {vec(bsxfun(@plus,  (0:1)', (0:7)*34)), ...
         vec(bsxfun(@plus,  (1:6)', (0:7)*34)), ...
         vec(bsxfun(@plus,  (8:13)', (0:7)*34)), ...
         vec(bsxfun(@plus,  (15:20)', (0:7)*34)), ...
         vec(bsxfun(@plus,  (22:27)', (0:7)*34)), ...
         vec(bsxfun(@plus,  (29:34)', (0:7)*34))};
idxTick = {[0, 0.5, 1, 1.5], [0,1,2], [0,2,4], [0,3,6,9], ...
           [0,5,10,15], [0,5,15,25]}; 

 [~, idxT] = min(abs(output{1}.Ts-1));

limitRate = zeros(length(idxSP),1); %
 % compute estimate of divergence rate 
dcN = zeros(length(idxSP),length(output{1}.Ns));    

for i = 2:length(idxSP)
if splitFigs
  figure(30+i)
else
  subplot(21,34,idxSP{i})
end    


 for j = find(output{i}.Ns<=320, 1, 'last'):-1:1   
  plot(output{i}.Ts, output{i}.cN(j,:), 'color', clrs(j,:), ...
      'linewidth', 1.5, 'markerSize', 1.5); hold on
  lgnd{find(output{i}.Ns<=320, 1, 'last')-j+1} = ['n = ', ...
      num2str(output{i}.Ns(j))];
  axis([min(output{i}.Ts),max(output{i}.Ts), 0, ...
      1.05*max(vec(output{i}.cN(1:find(output{i}.Ns<=320, 1, 'last'),:)))])

 end
  line([1,1], [0, 1.05*max(output{i}.cN(:))],'lineStyle', '--', ...
      'color', 'k', 'linewidth', axesThickness)
 
  box off
  set(gca, 'TickDir', 'out')
  set(gca, 'FontSize', fontSize)
  xlabel('temperature', 'FontName', fontName, ...
      'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
  set(gca, 'XTick', [1, 1.5, 2])
  if i == 2
   ylabel('specific heat', 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
  end
  set(gca, 'YTick', idxTick{i})
  text(1.3, 0.9*max(vec(output{i}.cN(1:find(output{i}.Ns<=320,1, ...
      'last'),:))), ['\rho = ', num2str(output{i}.rho)], ...
      'FontName', fontName, 'FontSize', fontSizeXlabel, ...
      'FontWeight', fontWeight )
  set(gca, 'Linewidth', axesThickness)
 
  dcN(i,:) =  output{i}.cN(:,idxT)./output{i}.Ns';
  limitRate(i) = ( output{i}.a*(output{i}.a+1)*psi(1,output{i}.a+1) + ...
                  output{i}.b*(output{i}.b+1)*psi(1,output{i}.b+1) )/ ...
                ((output{i}.a+output{i}.b)*(output{i}.a+output{i}.b+1))...
                - psi(1,output{i}.a+output{i}.b+1) ...
                + output{i}.a*output{i}.b*(psi(0,output{i}.a+1)- ...
                  psi(0,output{i}.b+1))^2/...
                ((output{i}.a+output{i}.b)^2*...
                 (output{i}.a+output{i}.b+1));
end

legend(lgnd);

%% subplot 'b.0)' : plot c(T=1) for different rho's

 clrs = [
254,227,145;   %
254,196,79;    % a few 
254,153,41;    % more and 
236,112,20;    % slightly 
204,76,2;      % different
153,52,4;      % colors
102,37,6]/255; %


lgnd = cell(length(output)-1,1);
for i = length(output):-1:2 % do not use rho = 0.001, it is just too small
    [~, idxT] = min(abs(output{i}.Ts-1));
    plot(output{i}.Ns, output{i}.cN(:,idxT), 's-', 'color',  clrs(i,:), ...
        'linewidth', 3, 'markerSize', 4)
    hold on
    lgnd{end-i+2} = ['\rho = ', num2str(output{i}.rho)];
    output{i}.cN(:,idxT)
end
legend(lgnd, 'Location', 'Northwest')
legend boxoff
box off
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', fontSize)
xlabel('population size', 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
set(gca, 'YTick', [0,10,20,30])
ylabel('specific heat at T =  1', 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
set(gca, 'XTick', output{1}.Ns)

%% subplot b) Summarize C(T=1) as function of N and estimated limit N->Inf
if splitFigs
  figure(37)
  subplot(1,2,1)
else
  subplot(21,34,vec(bsxfun(@plus, (1:8)', (11:20)*34)))
end    

 clrs = [
254,227,145;   %
254,196,79;    % a few 
254,153,41;    % more and 
236,112,20;    % slightly 
204,76,2;      % different
153,52,4;      % colors
102,37,6]/255; %

  for i = 1:length(idxSP) % got to plot phantom lines (just one point) to 
                         % get figure legend to appear in the correct order
  semilogx(output{1}.Ns(1), dcN(i,1), '-', 'color', clrs(end+1-i,:), ...
      'linewidth', 0.001); hold on;      
  end 
  
  for i = 2:length(idxSP)
  semilogx(output{1}.Ns, dcN(i,:), '-', 'color', clrs(i,:), ...
      'linewidth', 2); hold on;
  semilogx(output{1}.Ns(end)*[2,2.2], limitRate(i)*[1,1], ...
      'color', clrs(i,:), 'markerSize', 2, 'linewidth', 2);
  semilogx(output{1}.Ns, dcN(i,:), 's', 'color', clrs(i,:), ...
      'linewidth', 1.5, 'markerSize', 1.5); 
  end
 set(gca, 'TickDir', 'out'); box off
 set(gca, 'Linewidth', axesThickness)
 
 set(gca, 'FontSize', fontSize)
 xlabel('population size', 'FontName', fontName, ...
     'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
 ylabel('sp. heat divergence rate', 'FontName', fontName, ...
     'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
 ticks = [output{1}.Ns, output{1}.Ns(end)*2.05];
 set(gca,'XTick', ticks)
 lbls = cell(length(ticks),1);
 for i = 1:length(ticks), 
     lbls{i} = num2str(ticks(i)); 
 end; lbls{end} = '';
 lgnd = cell(length(idxSP),1);
 for i = 2:length(idxSP)
     lgnd{i} = ['\rho = ', num2str(output{end+1-i}.rho)];
 end
 if isempty(lgnd{1})
     lgnd = lgnd(2:end);
 end
 hold off
 legend(lgnd)
 set(gca,'XTickLabel', lbls)

 set(gca, 'YTick', 0.01:0.01:0.05);
 axis([20^0.9, output{i}.Ns(end)*2.2, 0.8*min(dcN(:)), 1.05*max(dcN(:))])

% some voodoo to replace the last tick with the symbol for infinity 
xTicks = get(gca, 'xtick'); yTicks = get(gca, 'ytick');
ax = axis; %Get left most x-position
HorizontalOffset = 0.1;
% Reset the xtick labels in desired font 
minY = min(yTicks);
text(output{i}.Ns(end)*2.05, 0.8*min(dcN(:)) - 0.0015, ...
    ' $ \mathbf{\infty}$',...
'HorizontalAlignment','Right','interpreter', 'latex');   


% subplot b) Summarize C(T=1)/n as function of rho
if splitFigs
  figure(37)
  subplot(1,2,2)
else
  subplot(21,34,vec(bsxfun(@plus, (9:11)', (11:20)*34)))
end    

rhos = zeros(length(output), 1);
for i = 2:length(rhos), rhos(i) = output{i}.rho; end;

plot(rhos(:), limitRate(:), '-', 'color', 0 * [1,1,1]), hold on
for i = 2:length(idxSP),
  plot(rhos(i), limitRate(i), 's', 'color', clrs(i,:), ...
      'markerSize', 1.5, 'linewidth', 1.5)
end

plot(rhos, rhos .* log( (1-output{1}.mu)/output{1}.mu )^2 * ...
    ( (1-output{1}.mu) * output{1}.mu ), 's', 'color', 0.3*[1,1,1], ...
    'linewidth', 1.0, 'markerSize', 1.5, 'linewidth', 1.5)
plot(rhos, rhos .* log( (1-output{1}.mu)/output{1}.mu )^2 * ...
    ( (1-output{1}.mu) * output{1}.mu ), '--', 'color', 0*[1,1,1], ...
    'linewidth', 1.0, 'linewidth', 1)

plot(rhos(:), limitRate(:), '-', 'color', 0 * [1,1,1]), hold on
for i = 2:length(idxSP),
  plot(rhos(i), limitRate(i), 's', 'color', clrs(i,:), ...
      'markerSize', 1.5, 'linewidth', 1.5)
end

hold off

set(gca, 'FontSize', fontSize)
xlabel('correlation', 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
box off
set(gca, 'XTick', [0, 0.1, 0.2])
set(gca, 'YTick', [0, 0.1, 0.2])
set(gca, 'TickDir', 'out')
axis([0, 1.05*max(rhos), 0.8*min(dcN(:)), 1.05*max(dcN(:))])

%% subplot c) Heat traces for different stimuli

maxi = 6;

clrs = [254,153,41;
236,112,20;
204, 76,2;      % colors for inidivual traces.
153,52,4;
102,37,6;
0,0,0]/255;


for j = 1:3
if j==1    
 if splitFigs
  figure(38)
  subplot(1,3,2)
 else
  subplot(21,34,vec(bsxfun(@plus, (21:26)', (10:20)*34)))
 end    
 load('../results/K_pairwise_final/heat_traces_nat.mat') % specific heat
 rho = 0.075;
 idxRepetMax = 10;
 yTick = [1,2,3,4]; %[0.4, 0.8, 1.2];
elseif j==2
if splitFigs
  figure(38)
  subplot(1,3,1)
else
  subplot(21,34,vec(bsxfun(@plus, (28:33)', (10:20)*34)))
end 
 rho = 0.341;
 load('../results/K_pairwise_final/heat_traces_fff.mat') % specific heat
 idxRepetMax = 10;
elseif j ==3
if splitFigs
  figure(38)
  subplot(1,3,3)
else
  subplot(21,34,vec(bsxfun(@plus, (14:19)', (10:20)*34)))
end
 rho = 0.033;
 load('../results/K_pairwise_final/heat_traces_cb.mat') % specific heat
 idxRepetMax = 10;
end
maxi = size(cN,1)/2;
lgnd = cell(maxi,1);
for i = maxi:-1:1, 
  plot(Ts, squeeze(mean(cN(2*i,:, :))), '-', ...
      'color', clrs(i,:), 'linewidth', 2.5); 
  lgnd{i} = ['n = ', num2str(20*(maxi-i+1))];
  hold on;  
  plot(Ts, squeeze(cN(2*i,:, :)), '-', ...
      'color', clrs(i,:), 'linewidth', 1); 
  lgnd{i} = ['n = ', num2str(20*(maxi-i+1))];
  
end;
line([1,1], [0, 4.5], 'lineStyle', '-', 'color', 'k', 'linewidth', 0.8)
for i = maxi:-1:1,
  plot(Ts, squeeze(mean(cN(2*i,:, :))), '-', ...
      'color', clrs(i,:), 'linewidth', 2.5); 
  hold on;  
  plot(Ts, squeeze(cN(2*i,:, :)), '-', ...
      'color', clrs(i,:), 'linewidth', 1); 
  lgnd{i} = ['n = ', num2str(20*(maxi-i+1))];
  
end;
for i = maxi:-1:1,
    plot(Ts, squeeze(mean(cN(2*i,:, :))), '-', ...
        'color', clrs(i,:), 'linewidth', 0.8);
    hold on;
  plot(Ts, squeeze(cN(2*i,:, :)), '-', ...
      'color', clrs(i,:), 'linewidth', 1); 
  lgnd{i} = ['n = ', num2str(20*(maxi-i+1))];
  
end
if j == 1
  legend(lgnd, 'location', 'East')
end
hold off

set(gca, 'Linewidth', axesThickness)
text(1.3, 0.9*4.5, ['\rho = ', num2str(rho)], ...
    'FontName', fontName, 'FontSize', fontSizeXlabel, ...
    'FontWeight', fontWeight )
box off
set(gca, 'TickDir', 'out')
set(gca, 'XTick', [1, 1.5, 2])
set(gca, 'FontSize', fontSize)
if j == 3
 set(gca, 'YTick', yTick)
 ylabel('specific heat c', 'FontName', fontName, ...
     'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
else
 set(gca, 'YTick', [])
end    
xlabel('temperature', 'FontName', fontName, ...
    'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
axis([0.8, 2, 0, 4.5])
end


% %% export figure to specified size
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
%  hgexport(figure3,'f1.pdf',s);     
%  export_fig('fig3.pdf')
%  
%  
%  
% %              s.Width= '3';
% %             s.Height= '2';
%              
% % hgexport(inset3_1,'test.pdf',s);     
% % export_fig('fig3_1.pdf')
%  
% end