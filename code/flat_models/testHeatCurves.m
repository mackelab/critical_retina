clc
clear all
%close all

 % Load P(K) from simulation data:
 load('C:\Users\Loki\Desktop\Criticality\code_critical_retina\Model Retina\Retina_model_data.mat')
 n = size(output,3);
 pcount = sum(output,3);
 pcount = histc(pcount(:),0:n);
 pcount = pcount/sum(pcount);
 % so, we actually only want to keep pcount to avoid clutter:
 clear cell_positions center_surround_ratio_large center_surround_ratio_small
 clear firing_rate_per_bin noise_radius noise_ratio pixel_size
 clear sigma_large_RF sigma_small_RF small_RF_mask stimulus
 clear threshold_offset threshold_slope trials weigth_surr
 clear white_noise_std output n
% Various other sourves of P(K):
%
% Salamander data N=40
%load('C:\Users\Loki\Desktop\Criticality\Code\Package\Data\SalamanderDataStatistics.mat');
%pcount = data.populationAcitivityCounts_valuesForKfrom0to40; % get true count distribution
% Data from tkacik 2014 'Searching for collective behavior' paper:
%load('C:\Users\Loki\Desktop\Criticality\Data\tkacik2014N100logPK.mat');
%pcount = zeros(101,1);
%pcount(1:46) = 10.^(tkacik2014N100logPK(:,2)');
%pcount = pcount/sum(pcount); % before this, sum(pcount) = 1.0169 
%load('C:\Users\Loki\Desktop\Criticality\Results\logGaussianPK100HeatCurves.mat');
% Data from previous run (to check for iterative criticality)
%load('C:\Users\Loki\Desktop\Criticality\Results\retSimN300_T1_0_HeatCurves.mat');
%[~, idxT]  = min(abs(Ts - 1));
%pcount = squeeze(pcountT(:, length(Ns), idxT )); clear idxT;
%pcount = 10.^(randn([101,1])-3); 
%pcount = pcount/sum(pcount);
%figure; plot(log(pcount)/log(10))
%pcount = binopdf(0:40, 40, 0.025)';
%pcount = nbinpdf(0:100, 1.5, 0.1)';
%pcount = pcount/sum(pcount);
%pcount = poisspdf(0:100, 10)';
%pcount = pcount/sum(pcount);
%pcount = 1./(1:335)';   % non-zero everywhere...
%pcount(floor(1*end/5):end) = 0; 
%pcount = pcount/sum(pcount);
%pcount = normpdf(linspace(0,1,321), 0.3, 0.1)';

%mu = 1.2/50; rho = 0.1;                     % Dichotomized Gaussian
%[~,~,pcount,~]=fit_flat_dg(mu,rho,333+1);
%pcount = pcount';


n = length(pcount)-1; % number of neurons

Ts = (0.9 : 0.00125 : 1.5)'; % range of temperatures

Ns = [20:20:320];         % range of network sizes (max(Ns) should be below 
                        % size(pcount) = n+1, of course)
                        
alphaOrT = 'T';     % decide whether to plainly vary the temperature T 
                    % or to do so while keeping the expected value of K
                    % constant by varying the 'alpha' in a related problem

clrs = hsv(length(Ns));       
idxRep = 10;
%--------------------------------------------------------------------------
pcountT = zeros(max(length(pcount)),length(Ns),length(Ts)); % P(K) at temperature T
varE = zeros(length(Ts),length(Ns)); % storages for variance of energy/T^2 
CNx = zeros(size(varE));             % and specific heat for each condition
CNk = zeros(size(varE));             % computed from x-space and K-space 
varLogPr = zeros(length(Ts),length(Ns)); %
varK = zeros(length(Ts),length(Ns));     % other stuff that might be interesting
h = zeros(length(Ts), length(Ns));       %

lgnd = cell(length(Ns),1);   % legend for figure

lambdaT1 = -Inf * ones(n+1,length(Ns)); % lambda at T=1 for various N
lambdaT1(1,:) = 0;

nSamples = 1000000; % may take a few seconds...
nHack = 1; % number of guaranteed samples for each K>0 to sneak into the data (makes sure P(K)>0)
s = sample_flat_model(pcount',nSamples-nHack*n); % draw full samples from count distribution
s(:,end+1:end+n*nHack) = repmat(triu(ones(n)),1,nHack);
s = sparse(s);

for j = 1:length(Ns)
  disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);
  hist_s = zeros(Ns(j)+1,1);              
  for t = 1:idxRep
   idxN = randsample(n,Ns(j));    % pick subset of Ns(j) <= n neurons
   histCounts = full(sum(s(idxN,:),1));          
   hist_s = hist_s + histc(histCounts,0:Ns(j))';   
  end
  hist_s = hist_s/sum(hist_s);
  
  %[~,~,hist_s,~]=fit_flat_dg(mu, rho, Ns(j));  % DG case without sampling issues
  %hist_s = hist_s';                            %
  
  %hist_s = binopdf( (0:Ns(j))', Ns(j), 0.025); % bernoulli case without sampling issues

  lambdaT1((1:Ns(j))+1,j) = fit_flat_maxent_model(hist_s); % get lambda from count distr.
 
  lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - gammaln(Ns(j)+1-(0:Ns(j))))';

 
 for i = 1:length(Ts)
   disp(['- temperature run #',num2str(i),' out of ',num2str(length(Ts))]);
   T = Ts(i);
   lambdaT = lambdaT1(:,j)/T;
  switch alphaOrT
   case 'T'
   logpcountT = -Inf*ones(n+1,1);
   logpcountT(1:Ns(j)+1) = lambdaT(1:Ns(j)+1) + lognchoosek;   
   pcountT(:,j,i) = exp(logpcountT);
   pcountT(:,j,i) = pcountT(:,j,i)/sum(pcountT(:,j,i)); % getting Z is easy here
   h(i,j) = 0;
   case 'alpha'
   [pcountT(1:(Ns(j)+1),j,i),~,models,h(i,j)]=flat_model_trace_temp(hist_s,1./T,true);
   varLogPr(i,j) = models.var_log_probs;
   pcountT(1:(Ns(j)+1),j,i) = pcountT(1:(Ns(j)+1),j,i) / sum(pcountT(1:(Ns(j)+1),j,i));
  end
   % Compute variance of E via sampling from K-space:
   tmp = cumsum([0;squeeze(pcountT(1:(Ns(j)+1),j,i))]); % works pretty quick
   [~,K] = histc(rand(nSamples,1), tmp);
   K = K-1;               % indexing...
   varK(i,j) = var(K);
  switch alphaOrT
   case 'T'
   E = lambdaT(K+1); % note that E already is lambdaT1 / T, i.e. var(E) = 1/T^2 * var(energy)      
   case 'alpha'
   E = (lambdaT(K+1) + h(i,j) * K); 
  end
   varE(i,j) = var(E); % var(energy)/T^2
   CNx(i,j) = (1/Ns(j)) * varE(i,j); % specific heat, up to Boltzmann constant k_B
   % Compute variance of log-probs directly from P(K) using Jakob's code (same results)
%   [varLogPr(i,j),~]= flat_model_var_log_probs(squeeze(pcountT(1:(Ns(j)+1),j,i)),exp(1));
%   CNk(i,j) = (1/Ns(j)) * varLogPr(i,j);
   
 end
 
end

clear E K CN T b histCounts hist_s i j t idx idxN lognchoosek logpcountT 

% save('C:\Users\Loki\Desktop\Criticality\Results\retSimN300_T1_0_alphaHeatCurves', ...
%'Ts', 'Ns', 'lambdaT1', 'n', 'nSamples', 'pcount', 'pcountT', 'varE',
%'varK', 'alphaOrT', 'h');

%% Preview figure for heat curves
textdx = .05; textdy = .1; % position of text relative to arrows (see fig)
for CNtype = 1
 figure;
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

%% Plot P( K | N, T, \lambda) for all used N, T (watch it become delta peaks)
CNtype = 1;
figure; 
for j = 1:length(Ns)
 for i = 1:length(Ts)
  tmp = squeeze(pcountT(:,j,i));
  plot(tmp); hold on;
  [m,v] = calc_mean_var(tmp);
  text(0.75*Ns(j), 0.5*max(tmp), ['var[x]: ', num2str(v)]);
  text(0.75*Ns(j), 0.3*max(tmp), ['E[x]: ',num2str(m)] ); hold off;
  title(['T = ', num2str(Ts(i)), ', N = ', num2str(Ns(j))])
  if strcmp(alphaOrT, 'alpha')
   lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - gammaln(Ns(j)+1-(0:Ns(j))))';  
   lambdaT = lambdaT1(:,j)/Ts(i);
   logpcountT = -Inf*ones(Ns(j)+1,1);
   logpcountT(1:Ns(j)+1) = lambdaT(1:(Ns(j)+1)) + h(i,j)*(0:Ns(j))' + lognchoosek;     
   tmp = exp(logpcountT); tmp = tmp/sum(tmp); 
   hold on; 
   plot(tmp, 'r--'); 
   hold off
  end
  pause
 end
end

%% Have a look at var(K) given N, T
figure(CNtype+45);
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
