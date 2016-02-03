% This script loads the spiking data of n simulated neurons, computes P(K) 
% from the spikes and then subsamples a new draw of patterns x (all neurons 
% are assumed identical w.r.t. firing rates and pairwise correlations) from
% P(K) to simulate activity at network sizes Ns(i)<=n, i = 1,...,length(Ns)

clc
clear all
%close all

n = 1280;

mu = 1.5/50; 

rhos = [0.005, 0.01, 0.025, 0.05, 0.1, 0.25];
output = cell(length(rhos), 1);


Ts = (0.8 : 0.00125 : 2)'; % range of temperatures to look at 

Ns = [20, 40, 80, 160, 320, 640, 1280];

alphaOrT = 'T';     % decide whether to plainly vary the temperature T 
                    % or to do so while keeping the expected value of K
                    % constant by varying the 'alpha' in a related problem.
                    % Leave at 'T' if unsure, or consult Tkacik et al. 2014
                    % for details of varying said 'alpha'

ifApproxByBinomial     = false;                    
ifApproxByBetaBinomial = true; 
ifApproxByFlatIsing    = false;

ifPlot = false;

nHack = 1; % number of guaranteed samples for each K>0 to sneak into 
           % the data (makes sure P(K)>0)

           
for r = 1:length(rhos);


          
rho = rhos(r);                     

% This is alternative code for computing the \alpha and \beta parameters 
% from beta-binomial fits to data drawn from a dichotomized Gausisan
% with specified mean mu and correlation rho. 
%        [~,~,pcount,~] = fit_flat_dg(mu,rho,Ns(end));
%        pcount = pcount';        
%            mu1n =  (0:Ns(end)) * pcount;
%            mu2n = (0:Ns(end)).^2 * pcount;
%            Z = ( Ns(end) * (mu2n/mu1n - mu1n -1)) + mu1n;
%            output{r}.a = (Ns(end) * mu1n - mu2n) / Z;
%            output{r}.b = (Ns(end) - mu1n) * (Ns(end) - mu2n/mu1n) / Z;
           
output{r}.a =    mu    * (1/rho -1);                        
output{r}.b =  (1 -mu) * (1/rho -1);

%--------------------------------------------------------------------------

pcountT = zeros(max(n+1),length(Ns),length(Ts)); % P(K) at temperature T
varE = zeros(length(Ns),length(Ts)); % storages for variance of energy/T^2 
cN = zeros(size(varE));             % and specific heat for each condition
varLogPr = zeros(length(Ts),length(Ns)); %
h = zeros(length(Ts), length(Ns));       %
a = zeros(1, length(Ns)); 
b = zeros(1, length(Ns)); 

%mu1 = zeros(length(Ns));
%mu2 = zeros(length(Ns));

explogpcount = zeros(Ns(end)+1, length(Ns));
pcounts = zeros(size(explogpcount));

lgnd = cell(length(Ns),1);   % legend for figure

lambdaT1 = -Inf * ones(n+1,length(Ns)); % lambda at T=1 for various N
lambdaT1(1,:) = 0;

for j = 1:length(Ns)
    disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);
                    
        if ifApproxByBinomial
            [~,~,hist_s,~] = fit_flat_dg(mu,rho,Ns(j));
            hist_s = hist_s';            
            hist_s = binopdf(0:Ns(j), Ns(j), mean(histCounts)/Ns(j));
        end
        
        if ifApproxByBetaBinomial

            a(j) = output{r}.a; % does not depend on the actual current
            b(j) = output{r}.b; % subpopulation size
            
% This is alternative code for computing the \alpha and \beta parameters 
% from beta-binomial fits to data drawn from a dichotomized Gausisan
% with specified mean mu and correlation rho. 
%            [~,~,hist_s,~] = fit_flat_dg(mu,rho,Ns(j));
%            hist_s = hist_s';            
%            mu1(j) =  (0:Ns(j)) * hist_s;
%            mu2(j) = (0:Ns(j)).^2 * hist_s;
%            Z = ( Ns(j) * (mu2(j)/mu1(j) - mu1(j) -1)) + mu1(j);
%            a(j) = (Ns(j) * mu1(j) - mu2(j)) / Z;
%            b(j) = (Ns(j) - mu1(j)) * (Ns(j) - mu2(j)/mu1(j)) / Z;
            lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - gammaln(Ns(j)+1-(0:Ns(j))))'; 
            logpcount = lognchoosek + betaln(a(j) + (0:Ns(j)), Ns(j) + b(j) - (0:Ns(j)))' - betaln(a(j), b(j));
            explogpcount(1:Ns(j)+1, j) = exp(logpcount);            
            hist_s = exp(logpcount);
        end
        
        
        if ifApproxByFlatIsing
            [~,~,hist_s,~] = fit_flat_dg(mu,rho,Ns(j));
            hist_s = hist_s';                        
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
                    varE (j,i) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp.^2) - (pcountT(1:(Ns(j))+1,j,i)' * tmp)^2;
                    cN(j,i) = (1/Ns(j)) * varE(j,i); % specific heat, up to Boltzmann constant k_B
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

clear E K T histCounts hist_s logpcount i j t idx lognchoosek logpcountT 


%% Preview figure for heat curves
if ifPlot
%figure(42)
subplot(1,6,(7:9)-6)
clrs = jet(length(Ns));       % colors used by preview plot
ct = 1;
 for j = 1:1:length(Ns)
   plot(Ts, squeeze(cN(j,:)), 'color', clrs(j,:), 'marker', 'x', 'linewidth', 1.5, 'markerSize', 6); hold on;
  lgnd{ct} = ['N = ', num2str(Ns(j))];
  ct = ct + 1;
 end
 axis([min(Ts)-0.1,max(Ts)+0.1,0.9*min(cN(:)),1.05*max(cN(:))])
 xlabel('temperature T'), ylabel('specific heat var[E]/(k_B T^2)')
 box off; set(gca, 'TickDir', 'out')
 legend(lgnd);
 hold off;
 title(['Heat curves for DG activity data with N = ', num2str(Ns(end)), ' cells, rho =', num2str(rho)])

 clear pieakHi tmp 

subplot(1,6,(10:12)-6)
[~, idxT] = min(abs(Ts - 1));
plot(Ns, squeeze(cN(:,idxT))', 'bo-', 'linewidth', 2)
hold on
plot(Ns, squeeze(max(cN(:,:), [],2))', 'ko-', 'linewidth', 2)
hold off
%line([Ns(1), Ns(end)], [cN(idxT, 1), cN(idxT, end)], 'color', 'b')     % reference lines (will be very similar
%line([Ns(1), Ns(end)], [max(cN(:, 1)), max(cN(:, end))], 'color', 'k') % to the results from random subsampling)
legend('c(T=1)', 'max\{c(T)\}')
title('specific heat divergence against system size n')
xlabel('system size n')
ylabel('c(T)')

end

%%
output{r}.mu  = mu;
output{r}.rho = rho;
output{r}.cN = cN;
output{r}.Ns = Ns;
output{r}.Ts = Ts;
output{r}.nHack = nHack;
output{r}.pcounts = pcounts;
output{r}.lambdaT1 = lambdaT1;

%output{r}.mu1  = mu1;
%output{r}.mu2  = mu2;

output{r}.as = a;
output{r}.bs = b;
%pause; 
end
%%
%save('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Figures/Fig3/fig3_data.mat', ...
%     'output')