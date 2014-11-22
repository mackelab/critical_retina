clc
clear all
%close all

% load('C:\Users\Loki\Desktop\Criticality\code_critical_retina\Model Retina\Retina_model_data.mat')
% n = size(output,3);
% pcount = sum(output,3);
% pcount = histc(pcount(:),0:n);
% pcount = pcount/sum(pcount);
% clear cell_positions center_surround_ratio_large center_surround_ratio_small
% clear firing_rate_per_bin noise_radius noise_ratio pixel_size
% clear sigma_large_RF sigma_small_RF small_RF_mask stimulus
% clear threshold_offset threshold_slope trials weigth_surr
% clear white_noise_std output n
%load('C:\Users\Loki\Desktop\Criticality\Code\Package\Data\SalamanderDataStatistics.mat');
%pcount = data.populationAcitivityCounts_valuesForKfrom0to40'; % get true count distribution

% Data from tkacik 2014 'Searching for collective behavior' paper:
%load('C:\Users\Loki\Desktop\Criticality\Data\tkacik2014N100logPK.mat');
%pcount = zeros(101,1);
%pcount(1:46) = 10.^(tkacik2014N100logPK(:,2)');
%pcount = pcount/sum(pcount); % before this, sum(pcount) = 1.0169 

%load('C:\Users\Loki\Desktop\Criticality\Results\logGaussianPK100HeatCurves.mat');

% Data from previous run (to check for iterative criticality)
%load('C:\Users\Loki\Desktop\Criticality\Results\tkacik2014N100_T1_0_HeatCurves.mat');
%[~, idxT]  = min(abs(Ts - 0.9));
%pcount = squeeze(pcountT(:, length(Ns), idxT )); clear idxT;

%pcount = 10.^(randn([101,1])-3); 
%pcount = pcount/sum(pcount);
%figure; plot(log(pcount)/log(10))

pcount = binopdf(0:334, 334, 0.025)';

%pcount = nbinpdf(0:100, 1.5, 0.1)';
%pcount = pcount/sum(pcount);

%pcount = poisspdf(0:100, 10)';
%pcount = pcount/sum(pcount);

clear CNk CNx Ns Ts pcountT

n = length(pcount)-1; % number of neurons

Ts = (0.5 : 0.01 : 2.0)'; % range of temperatures
Ns = 50:50:300;         % range of network sizes (max(Ns) should be below 
                        % size(pcount) = n+1, of course)
clrs = hsv(length(Ns));       
idxRep = 10;
%--------------------------------------------------------------------------
pcountT = zeros(max(length(pcount)),length(Ns),length(Ts)); % P(K) at temperature T
varE = zeros(length(Ts),length(Ns)); % storages for variance of energy and
CNx = zeros(size(varE));              % specific heat for each condition
CNk = zeros(size(varE));              
varLogPr = zeros(length(Ts),length(Ns));
varK = zeros(length(Ts),length(Ns)); 

lgnd = cell(length(Ns),1);   % legend for figure
textdx = .05; textdy = .1; % position of text relative to arrows (see fig)

lambdaT1 = -Inf * ones(n+1,length(Ns));
lambdaT1(1,:) = 0;

nSamples = 1000000;
%s = sample_flat_model(pcount',nSamples); % draw full samples from count distribution
%s = sparse(s);

for j = 1:length(Ns)
  disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);
%  hist_s = zeros(Ns(j)+1,1);              
%  for t = 1:idxRep
%   idxN = randsample(n,Ns(j));    % pick subset of Ns(j) <= n neurons
%   histCounts = full(sum(s(idxN,:),1));          
%   hist_s = hist_s + histc(histCounts,0:Ns(j))';   
%  end
%  hist_s = hist_s/sum(hist_s);
  hist_s = binopdf( (0:Ns(j))', Ns(j), 0.025);
  lambdaT1((1:Ns(j))+1,j) = fit_flat_maxent_model(hist_s); % get lambda from count distr.
  
  lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - gammaln(Ns(j)+1-(0:Ns(j))))';

 
 for i = 1:length(Ts)
   disp(['- temperature run #',num2str(i),' out of ',num2str(length(Ts))]);
   T = Ts(i);
   lambdaT = lambdaT1(:,j)/T;
   logpcountT = -Inf*ones(n+1,1);
   logpcountT(1:Ns(j)+1) = lambdaT(1:Ns(j)+1) + lognchoosek;   
   pcountT(:,j,i) = exp(logpcountT);
   pcountT(:,j,i) = pcountT(:,j,i)/sum(pcountT(:,j,i)); % getting Z is easy here

   %[pcountT(1:(Ns(j)+1),j,i),~,models]=flat_model_trace_temp(hist_s,1./T,1);
   %varLogPr(i,j) = models.var_log_probs;
   %pcountT(1:(Ns(j)+1),j,i) = pcountT(1:(Ns(j)+1),j,i) / sum(pcountT(1:(Ns(j)+1),j,i));
   
   % Compute variance of E via sampling from space of full patterns:
   tmp = cumsum([0;squeeze(pcountT(1:(Ns(j)+1),j,i))]); 
   [~,K] = histc(rand(nSamples,1), tmp);
   K = K-1;               varK(i,j) = var(K);
   E = lambdaT(K+1);      varE(i,j) = var(E);
   CNx(i,j) = (1/Ns(j)) * varE(i,j); % specific heat, up to Boltzmann constant k_B
   % Compute variance of E directly from P(K) (it's a flat model after all)
%   [varLogPr(i,j),~]= flat_model_var_log_probs(squeeze(pcountT(1:(Ns(j)+1),j,i)),exp(1));
%   CNk(i,j) = (1/Ns(j)) * varLogPr(i,j);
   
 end
 
end

clear E K CN T b h histCounts hist_s i j t idx idxN lognchoosek logpcountT 

for CNtype = 1
 figure(CNtype+42)
 switch CNtype
   case 1 % var(E) via full space of x's
     CN = bsxfun(@times, CNx, 1); % CNx
   case 2 % var(E) directly from P(K)
     CN = bsxfun(@times, CNk, 1); %CNk;
 end
 for j = 1:length(Ns)
  plot(Ts, CN(:,j), 'color', clrs(j,:), 'marker', 'o', 'linewidth', 1.5, 'markerSize', 1.5); hold on;
 end
 boxl = 0.75*max(Ts); boxr = max(Ts); boxu = 2/3*max(CN(:)); boxd = 1/6*max(CN(:));
 for j = 1:length(Ns)

 [peakHi,peakLoc] = max(CN(:,j));
  line([Ts(peakLoc),Ts(peakLoc)],[peakHi,peakHi+0.25],'color','g','linewidth',1.5);
  line([Ts(peakLoc),Ts(peakLoc)+0.02],[peakHi,peakHi+0.1],'color','g','linewidth',1.5);
  line([Ts(peakLoc),Ts(peakLoc)-0.02],[peakHi,peakHi+0.1],'color','g','linewidth',1.5);
  text(Ts(peakLoc)+textdx, peakHi+textdy, ...
       ['T = ', num2str(Ts(peakLoc)), ', C/N = ', num2str(peakHi), ', C/N^2 = ', num2str(peakHi/Ns(j))], ...
        'color', 'g');
  lgnd{j} = ['N = ', num2str(Ns(j))];
  
  plot(boxl+0.1*(boxr-boxl)+0.8*(Ns(j))/(max(Ns))*(boxr-boxl), boxd+0.1*(boxu-boxd)+0.8*peakHi/max(CN(:))*(boxu-boxd), 'g*');

 end
 text(boxl+0.49*(boxr-boxl), boxd-0.1, 'N');
 text(boxl+0.17*(boxr-boxl), boxu+0.1, 'Dependence of peak C/N(T) on N');
 text(boxl-0.1, boxd+0.5*(boxu-boxd), 'C/N');
 line([1,1],[0,1.05*max(varE(:))], 'color', 'r', 'linewidth', 1.5)
 line([boxl,boxl],[boxu,boxd], 'color', 'k');
 line([boxr,boxr],[boxu,boxd], 'color', 'k');
 line([boxl,boxr],[boxu,boxu], 'color', 'k');
 line([boxl,boxr],[boxd,boxd], 'color', 'k');
 line([boxl,boxr],[boxd,boxu], 'color', 'b');
 axis([min(Ts)-0.1,max(Ts)+0.1,0.9*min(CN(:)),1.05*max(CN(:)+textdy)])
 xlabel('temperature T'), ylabel('specific heat var[E]/(k_B T^2)')
 box off; set(gca, 'TickDir', 'out')
 legend(lgnd);
 hold off;
 title('Heat curves for flat maxEnt model for population activity count data with N = 100 (Tkacic et al, PlosCB, Jan 2014)')

end
clear pieakHi tmp CNtype

%%
CNtype = 1;
figure(CNtype+43); 
for j = 1:length(Ns)
 for i = 1:length(Ts)
  plot(squeeze(pcountT(:,j,i))); hold on;
  text(0.75*n, 0.1, num2str(var(squeeze(pcountT(:,j,i)))));
  text(0.75*n, 0.07, num2str(var(squeeze(pcountT(:,j,i)))/Ns(j)^2)); hold off;
  title(['T = ', num2str(Ts(i)), ', N = ', num2str(Ns(j))])
  pause
 end
end

%% 
figure(CNtype+44);
subplot(3,3,[1:2,4:5,7:8])
varKpK = zeros(length(Ns), length(Ts));
for j = 1:length(Ns)
 for i = 1:length(Ts)
    %varK(j,:) = var(squeeze(pcountT(:,j,:)))';
    tmp = squeeze(pcountT(:,j,i));
    varKpK(j,i) = sum( ((0:n).^2)' .* tmp ) - sum( (0:n)' .* tmp).^2;
 end
end
imagesc(varKpK); colorbar
set(gca, 'YTickLabel', Ns);
set(gca, 'XTick', 1:10:length(Ts))
set(gca, 'XTickLabel', Ts(1:10:end));
xlabel('Temperature T'), ylabel('Network size N')
title('var(K) as a function of temperature and network size')
subplot(3,3,3)
[~,idx] = min(abs(Ts-0.8));
plot(Ns, varKpK(:,idx)); hold on
plot(Ns, varK(idx,:), 'k-')
b = glmfit([Ns',Ns'.^2], varKpK(:,idx), 'normal', 'constant', 'off');
plot(Ns,[Ns;Ns.^2]'*b, 'r'); 
aPar = b(1); bPar = -b(2)/(2*aPar); 
text(0.7*max(Ns), 0.3*max(varKpK(:,idx)), ['a =', num2str(aPar)])
text(0.7*max(Ns), 0.2*max(varKpK(:,idx)), ['b =', num2str(bPar)])
clear aPar bPar, hold off;
xlabel('Network size N'), ylabel('var(K) at T=1')
legend('var_N(K)', 'var_N(K) sampling', 'a*(N-b)^2', 'location', 'NorthWest')
title('var(K) as a function of N at T = 0.8')
box off; set(gca, 'TickDir', 'out'); legend boxoff
axis([Ns(1), Ns(end), 0.95*min(varKpK(:,idx)), 1.05*max(varKpK(:,idx))])
%set(gca, 'XTick', 2:2:length(Ns))
%set(gca, 'XTickLabel', Ns(2:2:end));
subplot(3,3,6)
[~,idx] = min(abs(Ts-1));
plot(Ns, varKpK(:,idx)); hold on
plot(Ns, varK(idx,:), 'k-')
b = glmfit([Ns',Ns'.^2], varKpK(:,idx), 'normal', 'constant', 'off');
plot(Ns,[Ns;Ns.^2]'*b, 'r'); 
aPar = b(1); bPar = -b(2)/(2*aPar); 
text(0.7*max(Ns), 0.3*max(varKpK(:,idx)), ['a =', num2str(aPar)])
text(0.7*max(Ns), 0.2*max(varKpK(:,idx)), ['b =', num2str(bPar)])
clear aPar bPar, hold off; 
xlabel('Network size N'), ylabel('var(K) at T=1')
legend('var_N(K)', 'var_N(K) sampling','a*(N-b)^2', 'location', 'NorthWest')
title('var(K) as a function of N at T = 1')
box off; set(gca, 'TickDir', 'out'); legend boxoff
axis([Ns(1), Ns(end), 0.95*min(varKpK(:,idx)), 1.05*max(varKpK(:,idx))])
%set(gca, 'XTick', 2:2:length(Ns))
%set(gca, 'XTickLabel', Ns(2:2:end));
subplot(3,3,9)
[~,idx] = min(abs(Ts-1.1));
plot(Ns, varKpK(:,idx)); hold on
plot(Ns, varK(idx,:), 'k-')
b = glmfit([Ns',Ns'.^2], varKpK(:,idx), 'normal','constant', 'off');
plot(Ns,[Ns;Ns.^2]'*b, 'r')
aPar = b(1); bPar = -b(2)/(2*aPar); 
text(0.7*max(Ns), 0.3*max(varKpK(:,idx)), ['a =', num2str(aPar)])
text(0.7*max(Ns), 0.2*max(varKpK(:,idx)), ['b =', num2str(bPar)])
clear aPar bPar, hold off; 
xlabel('Network size N'), ylabel('var(K) at T=1')
legend('var_N(K)', 'var_N(K) sampling','a*(N-b)^2', 'location', 'NorthWest')
title('var(K) as a function of N at T = 1.1')
box off; set(gca, 'TickDir', 'out'); legend boxoff
axis([Ns(1), Ns(end), 0.95*min(varKpK(:,idx)), 1.05*max(varKpK(:,idx))])
%set(gca, 'XTick', 2:2:length(Ns))
%set(gca, 'XTickLabel', Ns(2:2:end));

%%
% save('C:\Users\Loki\Desktop\Criticality\Results\binomialP_0_025_HeatCurves.mat', ...
%      'Ts', 'Ns', 'lambdaT1', 'n', 'nSamples', 'pcount', 'pcountT', 'varE', 'varK');