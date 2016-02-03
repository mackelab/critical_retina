%% produces figure four for the journal version of our 'criticality and correlations' project
% figure '4' is the figure summarising effects of subsampling schemes. 

clear all

addpath(genpath('../code/'))

N = 100;        % overall population size
ns = 5:5:N; % sizes of subpopulations that are explored
nsShow = [10:20:ns(end)]; % subpopulation sizes for which we plot c(T)
idxRep = 20;    % number of subpops for each size   

idxNr = cell(length(ns),1); % index of subpops ('r' for random, 'l' for 
idxNl = cell(length(ns),1); % linear sampling), should be fixed somewhere
for n = 1:length(ns)
    idxNr{n} = zeros(ns(n), idxRep);
    idxNl{n} = zeros(ns(n), idxRep);
    for i = 1:idxRep
        idxNr{n}(:,i) = randsample(N, ns(n));
        idxNl{n}(:,i) = 1:ns(n);
    end                                               
end

numSamples = 1000000; % number of samples drawn from dichotomized Gaussian
                      % (to add some noise to the correlations)

xDist = linspace(0,1,101);
Ts = (0.8 : 0.00125 : 2)'; % temperature range for specific heat curves
[~,idxT] = min(abs(Ts-1));


axesThickness  = 1.1; % overall thickness of axes. Doesn't seem to do much.
fontName = 'Arial';    fontWeight     = 'normal';
fontSize       = 1 * 10;   fontSizeTitle  = 1 * 16;   
fontSizeXlabel = 1 * 10;   fontSizeYlabel = 1 * 11;
fontSizeText   = 1 * 10;   fontSizeLegend = 1 * 11;
cmin = Inf; cmax = -Inf;
%% flat data case

% compute data for flat data case:
rho = 0.05;
% compute data for realistic data case:
Sigma = rho * ones(N,N);                  % set Sigma to desired values
Sigma(logical(diag(ones(N,1)))) = 1;      %
mu = 0.05*ones(N,1);
Cov = zeros(size(Sigma));                 % compute covariance matrix
idxRnd = randperm(N);
for i = 1:N
    for j = i:N
        Cov(i,j) = Sigma(idxRnd(i),idxRnd(j)) * (mu(1) * (1-mu(1)));
    end
end
Cov = Cov + Cov' - diag(diag(Cov));
X = sampleDichGauss01(mu,Cov,numSamples/100)'; % draw small sample
Sigma  = corrcoef(X');           % re-set Sigma to noisy empirical values
                                      
Cov = zeros(size(Sigma));                % re-compute covariance matrix
idxRnd = randperm(N); 
for i = 1:N
    for j = i:N
        Cov(i,j) = Sigma(idxRnd(i),idxRnd(j)) * (mu(1) * (1-mu(1)));
    end
end
Cov = Cov + Cov' - diag(diag(Cov));
X = sampleDichGauss01(mu,Cov,numSamples)'; % draw large sample
SigmaD = corrcoef(X');% 'empirical' correlation matrix
disp(['avg corr in flat data case is ', ...
    num2str(mean(SigmaD(~logical(diag(ones(N,1))))))])
cmin = min([cmin, min(SigmaD(~logical(diag(ones(N,1)))))]);
cmax = max([cmax, max(SigmaD(~logical(diag(ones(N,1)))))]);

%%
fig4=figure('Tag', 'fig4','units','normalized','position',[0,0,0.99,0.99]);
clrs = copper(length(ns));
brownIdx = 13;
% a)
subplot(3,4,1)
plot(xDist,rho*ones(size(xDist)),'color',clrs(brownIdx,:),'linewidth',2.5)
xlabel('RF centre distance')
ylabel('correlation')
box off
set(gca, 'TickDir', 'out')
axis([0,1,0,0.25])
set(gca, 'XTick', [0,0.5,1])
set(gca, 'YTick', [0,0.1,0.2])
set(gca, 'Linewidth', axesThickness)
% b) 
subplot(3,4,2)
imagesc(SigmaD-diag(diag(SigmaD)))
set(gca, 'XTick', [1,50,100])
set(gca, 'YTick', [1,50,100])
xlabel('neuron id')
ylabel('neuron id')
box off
set(gca, 'TickDir', 'out')
axis square
set(gca, 'Linewidth', axesThickness)
caxis([cmin, cmax])
% c) 
subplot(3,4,3)
tracen = zeros(length(ns),2);
  stdn = zeros(length(ns),2);
for n = 1:length(ns)
    tmp = zeros(idxRep,1);
    for i = 1:idxRep
        M = SigmaD(idxNr{n}(:,i),idxNr{n}(:,i));
        tmp(i) = mean(M(~logical(diag(ones(ns(n),1)))));
    end
    tracen(n,1) = mean(tmp);
    stdn(n,1)   = std(tmp);
    for i = 1:idxRep
        M = SigmaD(idxNl{n}(:,i),idxNl{n}(:,i));
        tmp(i) = mean(M(~logical(diag(ones(ns(n),1)))));
    end
    tracen(n,2) = mean(tmp);
    stdn(n,2)   = std(tmp);
end
plot(-10, -10, 'color', 0.4*[1,1,1],'linewidth', 2)
hold on
plot(-10, -10, 'color', clrs(brownIdx,:), 'linewidth', 2)
h = area(ns, [tracen(:,1)-stdn(:,1), 2*stdn(:,1)]); 
h(2).FaceColor = 0.6*[1,1,1]; h(1).FaceColor = [1,1,1];
h(2).EdgeColor = [1,1,1]; h(1).EdgeColor = [1,1,1];
plot(ns, tracen(:,1), 'color', 0.4*[1,1,1],'linewidth', 2)
plot(ns, tracen(:,2), 'color', clrs(brownIdx,:), 'linewidth', 2)
hold off
xlabel('population size')
ylabel('correlation')
axis([ns(1)-1,ns(end)+1,0.0,0.25])
set(gca, 'YTick', [0,0.1,0.2])
set(gca, 'XTick', [1,50,100])
box off
set(gca, 'TickDir', 'out')
set(gca, 'Linewidth', axesThickness)
% d) 
 cT1r = zeros(length(ns), idxRep);
 cTMr = zeros(length(ns), idxRep);
 cT1l = zeros(length(ns), idxRep);
 cTMl = zeros(length(ns), idxRep);
 for n = 1:length(ns)
     for i = 1:idxRep
         pcount = sum(X(idxNr{n}(:,i),:))';
         pcount = histc(pcount, 0:ns(n));
         pcount = pcount/sum(pcount);
         [cN, ~] = computeSpecificHeatTraces(pcount, Ts, 1);
 %        end
         cT1r(n,i) = cN(idxT);
         cTMr(n,i) = max(cN);
         pcount = sum(X(idxNl{n}(:,i),:))';
         pcount = histc(pcount, 0:ns(n));
         pcount = pcount/sum(pcount);
         [cN, ~] = computeSpecificHeatTraces(pcount, Ts, 1);
 %        end
         cT1l(n,i) = cN(idxT);
         cTMl(n,i) = max(cN);
     end
 end
% e) (inset!)
subplot(3,4,4)
plot(ns, mean(cT1r,2), 'color', 0.4*[1,1,1], 'linewidth', 2)
hold on
plot(ns, mean(cT1l,2), 'color', clrs(brownIdx,:), 'linewidth', 2)
hold off
xlabel('population size')
ylabel('specific heat')
axis([ns(1)-1, ns(end)+1, 0, 2.5])
box off
set(gca, 'XTick', [5,50,100])
set(gca, 'YTick', [0, 1, 2])
set(gca, 'TickDir', 'out')
set(gca, 'Linewidth', axesThickness)

%% spatial drop-off data case
clear X
% as this is just an illustration, we try to match the average corelation 
% as closely as possible - by hand. 
sigma2 = 1/8;    % given above numSamples, 
scale  = 0.2527; % works out empirically to give rho = 0.05 
f = @(x) scale * exp(-(1/sigma2* x).^2);
% compute data for realistic data case:
Sigma = f(ones(N,1)*(1:N)/N - (ones(N,1)/N*(1:N))');
Sigma(logical(diag(ones(N,1)))) = 1;
Cov = zeros(size(Sigma));
for i = 1:N
    for j = i:N
        Cov(i,j) = Sigma(i,j) * (mu(1) * (1-mu(1)));
    end
end
Cov = Cov + Cov' - diag(diag(Cov));
X = sampleDichGauss01(mu,Cov,numSamples)';
SigmaD = corrcoef(X');% 'empirical' correlation matrix
cmin = min([cmin, min(SigmaD(~logical(diag(ones(N,1)))))]);
cmax = max([cmax, max(SigmaD(~logical(diag(ones(N,1)))))]);
disp(['avg corr in realistic data case is ', ...
    num2str(mean(SigmaD(~logical(diag(ones(N,1))))))])

%%
% a)
subplot(3,4,5)
plot(xDist, f(xDist), 'linewidth', 2.5, 'color', clrs(brownIdx,:))
xlabel('RF centre distance')
ylabel('correlation')
box off
set(gca, 'TickDir', 'out')
axis([0,1,0,0.25])
set(gca, 'XTick', [0,0.5,1])
set(gca, 'YTick', [0,0.1,0.2])
set(gca, 'Linewidth', axesThickness)
% b) 
subplot(3,4,6)
imagesc(SigmaD-diag(diag(SigmaD)))
set(gca, 'XTick', [1,50,100])
set(gca, 'YTick', [1,50,100])
xlabel('neuron id')
ylabel('neuron id')
box off
set(gca, 'TickDir', 'out')
axis square
set(gca, 'Linewidth', axesThickness)
caxis([cmin, cmax])
subplot(3,4,2)
caxis([cmin, cmax])

% c) 
subplot(3,4,7)
tracen = zeros(length(ns),2);
  stdn = zeros(length(ns),2);
for n = 1:length(ns)
    tmp = zeros(idxRep,1);
    for i = 1:idxRep
        M = SigmaD(idxNr{n}(:,i),idxNr{n}(:,i));
        tmp(i) = mean(M(~logical(diag(ones(ns(n),1)))));
    end
    tracen(n,1) = mean(tmp);
    stdn(n,1)   = std(tmp);
    for i = 1:idxRep
        M = SigmaD(idxNl{n}(:,i),idxNl{n}(:,i));
        tmp(i) = mean(M(~logical(diag(ones(ns(n),1)))));
    end
    tracen(n,2) = mean(tmp);
    stdn(n,2)   = std(tmp);
end
plot(-10, -10, 'color', 0.4*[1,1,1],'linewidth', 2)
hold on
plot(-10, -10, 'color', clrs(brownIdx,:), 'linewidth', 2)
h = area(ns, [tracen(:,1)-stdn(:,1), 2*stdn(:,1)]); 
h(2).FaceColor = 0.6*[1,1,1]; h(1).FaceColor = [1,1,1];
h(2).EdgeColor = [1,1,1]; h(1).EdgeColor = [1,1,1];
plot(ns, tracen(:,1), 'color', 0.4*[1,1,1],'linewidth', 2)
plot(ns, tracen(:,2), 'color', clrs(brownIdx,:), 'linewidth', 2)
hold off
xlabel('population size')
ylabel('correlation')
axis([ns(1)-1,ns(end)+1,0.0,0.25])
set(gca, 'YTick', [0,0.1,0.2])
set(gca, 'XTick', [5,50,100])
box off
set(gca, 'TickDir', 'out')
set(gca, 'Linewidth', axesThickness)
% d) 
 cT1r = zeros(length(ns), idxRep);
 cTMr = zeros(length(ns), idxRep);
 cT1l = zeros(length(ns), idxRep);
 cTMl = zeros(length(ns), idxRep);
 for n = 1:length(ns)
     for i = 1:idxRep
         pcount = sum(X(idxNr{n}(:,i),:))';
         pcount = histc(pcount, 0:ns(n));
         pcount = pcount/sum(pcount);
         [cN, ~] = computeSpecificHeatTraces(pcount, Ts, 1);
         cT1r(n,i) = cN(idxT);
         cTMr(n,i) = max(cN);
         pcount = sum(X(idxNl{n}(:,i),:))';
         pcount = histc(pcount, 0:ns(n));
         pcount = pcount/sum(pcount);
         [cN, ~] = computeSpecificHeatTraces(pcount, Ts, 1);
         cT1l(n,i) = cN(idxT);
         cTMl(n,i) = max(cN);
     end
 end
% e) (inset!)
subplot(3,4,8)
plot(ns, mean(cT1r,2), 'color', 0.4*[1,1,1], 'linewidth', 2)
hold on
plot(ns, mean(cT1l,2), 'color', clrs(brownIdx,:), 'linewidth', 2)
hold off
xlabel('population size')
ylabel('specific heat')
axis([ns(1)-1, ns(end)+1, 0, 2.5])
set(gca, 'XTick', [5,50,100])
set(gca, 'YTick', [0, 1, 2])
box off
set(gca, 'TickDir', 'out')
set(gca, 'Linewidth', axesThickness)

%% add real data information

load('fig_data/fig5_data.mat') % output, idxC, pixel_size, cell_positions
  
clrs = hsv(3);

% a) data traces for spatial drop-off
subplot(3,4,9) 
for r = 1:3
    n = size(Output{r}.spikes,1);
    idxC = Res{r}.idxC;
    corrsSorted = Output{r}.spkCorrs(idxC(1:n),idxC(1:n));
    corrsSorted = corrsSorted - diag(diag(corrsSorted));
    Y = corrsSorted(logical(triu(ones(size(corrsSorted)),1)));
    dists = zeros(n);
    for i = 1:n
     for j = i+1 : n
      dists(i,j) = sum( (cell_positions(:,idxC(i))- ...
                         cell_positions(:,idxC(j)) ).^2 );
     end
    end
    dists = sqrt( dists );
    X = dists(logical(triu(ones(size(corrsSorted)),1)));
    x = floor(min(X)) : ceil(max(X));
    [~, idx] = histc(X, x);
    m = zeros(size(x)); v = zeros(size(x));
    for i = 1:length(x)
      tmp = Y(idx == x(i));
      m(i) = mean(tmp);
      v(i) = var(tmp);
    end
    plot(x, m, 'color', clrs(r,:), 'linewidth', 2)
    hold on;
end
legend('checkerboard', 'natural images', 'full-field flicker')
legend boxoff
set(gca, 'FontSize', fontSize)
set(gca, 'XTick', [0, 50, 100, 150])
set(gca, 'Linewidth', axesThickness)
xlabel('RF centre distance', 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
ylabel('correlation', 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
box off; set(gca, 'TickDir', 'out')
axis([0, 160, -0.02, 0.6])
set(gca, 'Ytick', [0,0.2,0.4])
hold off

for r = 1:3
    subplot(3,4,9+r)
    Ns = zeros(length(Res{r}.idxN),1);
    [~,idxT] = min(abs(Res{r}.Ts - 1 )); 
    for i = 1:length(Ns)
        Ns(i) = size(Res{r}.idxN{i},1);
    end
    for t = 1:size(Res{r}.idxN{i},2)
        tmp = squeeze(Res{r}.cNrnd(:,t,idxT))';
        plot(Ns, tmp, 'color', 0.4*[1,1,1], 'linewidth', 0.5)
        hold on
    end
    plot(Ns, squeeze(Res{r}.cNlin(:,1,idxT)), 'color', clrs(r,:), ...
        'linewidth', 2)   
    hold off
    set(gca, 'Linewidth', axesThickness)
    set(gca, 'FontSize', fontSize)
    set(gca, 'XTick', [0, 100, 200, 300])
    xlabel('population size', 'FontName', fontName, ...
        'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
    ylabel('specific heat at T = 1', 'FontName', fontName, ...
        'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
    box off; set(gca, 'TickDir', 'out')
    if r == 1
        set(gca, 'Ytick', [0,1,2])
        axis([16, 320, 0, 3])
        set(gca, 'XTick', [1,100,200,300])
    elseif r == 2
        set(gca, 'Ytick', [0,2,4])
        axis([16, 320, 0, 6])
        set(gca, 'XTick', [1,100,200,300])
    elseif r == 3
        set(gca, 'Ytick', [0,5,10])
        axis([16, 320, 0, 11])
        set(gca, 'XTick', [1,100,200,300])
    end        
    hold off
    
end


% %% export figure to specified size
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
%  hgexport(fig4,'f4.pdf',s);     
%  export_fig('fig4.pdf')
