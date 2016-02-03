% This script loads the spiking data of n simulated neurons, computes P(K) 
% from the spikes and then subsamples a new draw of patterns x (all neurons 
% are assumed identical w.r.t. firing rates and pairwise correlations) from
% P(K) to simulate activity at network sizes Ns(i)<=n, i = 1,...,length(Ns)

clc
clear all
%close all


load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_data_new2.mat')
load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Heat curves/idxSubsamples.mat')
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
clear i tmptmp tmp




n = 300;
idxRep = 10; 

Ts = (0.8 : 0.00125 : 2)'; % range of temperatures to look at 

Ns = [20:20:n];

alphaOrT = 'T';     % decide whether to plainly vary the temperature T 
                    % or to do so while keeping the expected value of K
                    % constant by varying the 'alpha' in a related problem.
                    % Leave at 'T' if unsure, or consult Tkacik et al. 2014
                    % for details of varying said 'alpha'

ifApproxByBetaBinomial = false; 

ifPlot = true;

nHack = 1; % number of guaranteed samples for each K>0 to sneak into 
           % the data (makes sure P(K)>0)

varE = zeros(length(Ns),length(Ts),idxRep,2); % storages for variance of energy/T^2 
cN = zeros(size(varE));             % and specific heat for each condition
           
for r = 1:2

           
 if r == 1
  ifApproxByBinomial     = true;                    
  ifApproxByFlatIsing    = false;
 elseif r == 2
  ifApproxByBinomial     = false;                    
  ifApproxByFlatIsing    = true;
 end
%--------------------------------------------------------------------------




pcountT = zeros(max(n+1),length(Ns),length(Ts)); % P(K) at temperature T
varLogPr = zeros(length(Ts),length(Ns)); %
%varK = zeros(length(Ts),length(Ns));     % other stuff that might be interesting
h = zeros(length(Ts), length(Ns));       %
a = zeros(1, length(Ns)); 
b = zeros(1, length(Ns)); 

explogpcount = zeros(Ns(end)+1, length(Ns));
pcounts      = zeros(Ns(end)+1, length(Ns));

lgnd = cell(length(Ns),1);   % legend for figure

lambdaT1 = -Inf * ones(n+1,length(Ns)); % lambda at T=1 for various N
lambdaT1(1,:) = 0;

for j = 1:length(Ns)
    disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);

    lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - gammaln(Ns(j)+1-(0:Ns(j))))';
        
  for t = 1:idxRep
        histCounts = full(sum(output.spikes(idxSubsamples{Ns(j)/10}(:,t),:),1));
        hist_s = histc(histCounts,0:Ns(j))';
        hist_s = hist_s/sum(hist_s);     
        
        if ifApproxByBinomial
            hist_s = binopdf(0:Ns(j), Ns(j), mean(histCounts)/Ns(j));
        end
        
        
        if ifApproxByFlatIsing
            mu1 =  (0:Ns(j)) * hist_s;
            mu2 = (0:Ns(j)).^2 * hist_s;
            [mu,rho,~]=meanvar_count_2_meancorr(mu1,mu2 - mu1^2,Ns(j));            
            [~,~,hist_s,~]=fit_flat_ising_model(mu,rho,Ns(j));
%            [hist_s,~]=flat_ising_count_distrib(hFI,2*JFI,Ns(j));
        end
        
        pcounts(1:Ns(j)+1, j) = hist_s; % store
        
        
        lambdaT1((1:Ns(j))+1,j) = fit_flat_maxent_model(hist_s); % get lambda from count distr.
        
        
        
        for i = 1:length(Ts)
            T = Ts(i);
            lambdaT = lambdaT1(:,j)/T;
            
            switch alphaOrT
                case 'T'
                    logpcountT = -Inf*ones(n+1,1);
                    logpcountT(1:Ns(j)+1) = lambdaT(1:Ns(j)+1) + lognchoosek;
                    pcountT(:,j,i) = exp(logpcountT);
                    pcountT(:,j,i) = pcountT(:,j,i)/sum(pcountT(:,j,i)); % getting Z is easy here
                    h(i,j) = 0;
                    tmp = lambdaT(1:(Ns(j)+1));
                    tmp(pcountT(1:(Ns(j))+1,j,i)==0) = 0;
                    varE (j,i,t) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp.^2) - (pcountT(1:(Ns(j))+1,j,i)' * tmp)^2;
                    cN(j,i,t,r) = (1/Ns(j)) * varE(j,i,t); % specific heat, up to Boltzmann constant k_B
                case 'alpha'
                    [pcountT(1:(Ns(j)+1),j,i),~,models,h(i,j)]=flat_model_trace_temp(hist_s,1./T,true);
                    varLogPr(i,j) = models.var_log_probs;
                    pcountT(1:(Ns(j)+1),j,i) = pcountT(1:(Ns(j)+1),j,i) / sum(pcountT(1:(Ns(j)+1),j,i));
                    % Compute variance of E via sampling from K-space:
                    tmp = cumsum([0;squeeze(pcountT(1:(Ns(j)+1),j,i))]); % works pretty quick
                    [~,K] = histc(rand(nSamples,1), tmp);
                    K = K-1;               % indexing...
                    varK(i,j) = var(K);
                    E = (lambdaT(K+1) + h(i,j) * K);
                    varE(i,j,t) = var(E,t); % var(energy)/T^2
                    cN(i,j,t) = (1/Ns(j)) * varE(i,j,t); % specific heat, up to Boltzmann constant k_B
            end
        end    
  end
end

clear E K T b histCounts hist_s logpcount i j t idx lognchoosek logpcountT 


%% Preview figure for heat curves
if ifPlot

subplot(2,2,r)
clrs = jet(length(Ns));       % colors used by preview plot
ct = 1;
 for j = 1:length(Ns)
  for t = 1:idxRep
   plot(Ts, squeeze(cN(j,:,t,r)), 'color', clrs(j,:), 'marker', 'x', 'linewidth', 1.5, 'markerSize', 6); hold on;
  end
%  lgnd{ct} = ['N = ', num2str(Ns(j))];
  ct = ct + 1;
 end
 axis([min(Ts)-0.1,max(Ts)+0.1,0.9*min(cN(:)),1.05*max(cN(:))])
 xlabel('temperature T'), ylabel('specific heat var[E]/(k_B T^2)')
 box off; set(gca, 'TickDir', 'out')
% legend(lgnd);
 hold off;
 title(['Heat curves for flat maxEnt model for simulated activity count data with N = ', num2str(Ns(end)), ' cells'])

 clear pieakHi tmp 

subplot(2,2,2+r)
[~, idxT] = min(abs(Ts - 1));
plot(Ns, squeeze(cN(:,idxT,:,r))', 'bo-', 'linewidth', 2)
hold on
for t = 1:idxRep
 plot(Ns, squeeze(max(squeeze(cN(:,:,t,r)), [],2))', 'ko-', 'linewidth', 2)
end
%line([Ns(1), Ns(end)], [cN(idxT, 1), cN(idxT, end)], 'color', 'b')     % reference lines (will be very similar
%line([Ns(1), Ns(end)], [max(cN(:, 1)), max(cN(:, end))], 'color', 'k') % to the results from random subsampling)
legend('c(T=1)', 'max\{c(T)\}')
title('specific heat divergence against system size n (should saturate)')
xlabel('system size n')
ylabel('c(T)')

end

%%

%pause; 
end
%save('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Figures/Fig2/fig2_data.mat', ...
%     'Ts', 'Ns', 'cN')