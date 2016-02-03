% This script loads the spiking data of n simulated neurons, computes P(K) 
% from the spikes and then subsamples a new draw of patterns x (all neurons 
% are assumed identical w.r.t. firing rates and pairwise correlations) from
% P(K) to simulate activity at network sizes Ns(i)<=n, i = 1,...,length(Ns)

clc
clear all
%close all

% Load P(K) from simulation data: 
% The location on your hard-drive has to be changed accordingly!
%load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_gaussian_neu_offset_029.mat')
load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_data_new2.mat')
load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Heat curves/idxSubsamples.mat')
disp('getting spatial subsampling order (may take a few seconds)')

%idxC = find(cell_positions(2, :) < 0);
idxC = 1:size(cell_positions,2);
cell_positions = cell_positions(:, idxC); 
output = output(:,:,idxC); 

getSpatialSubsamplingOrder; % cleans up unused variables, 
                            % gives re-ordered indices idxC,
                            % also re-arranges output into struct variable,
                            % and starts plotting results

%output.spikes = (rand([size(output.spikes,1), 1000000]) < 0.03);

n = size(cell_positions,2); % number of neurons in this simulation 

Ts = (0.7 : 0.00125 : 1.2)'; % range of temperatures to look at 

n = 316; 
%Ns = [15,30,60,120,240];% range of network sizes to look at (max(Ns) should 
Ns = [20:20:300];          % be below size(pcount) = n+1, of course)
                        
alphaOrT = 'T';     % decide whether to plainly vary the temperature T 
                    % or to do so while keeping the expected value of K
                    % constant by varying the 'alpha' in a related problem.
                    % Leave at 'T' if unsure, or consult Tkacik et al. 2014
                    % for details of varying said 'alpha'

ifApproxByBinomial     = false;                    
ifApproxByBetaBinomial = false; 
ifApproxByFlatIsing    = false;

ifRandomSubsampling    = true; % replace 'spatially informed' with random
idxRep = 50; % number of subsets drawn for each n to estimate P(K) at n

nHack = 1; % number of guaranteed samples for each K>0 to sneak into 
           % the data (makes sure P(K)>0)

%--------------------------------------------------------------------------
if  ~ifRandomSubsampling 
 idxRep = 1; % does not make much sense to have this > 1 in this case
end

% From here on, the proper algorithm starts (no need to change stuff)
 pcount = sum(output.spikes,1)';
 pcount = histc(pcount,0:n);
 pcount = pcount/sum(pcount);

 idxShuffle = randsample(n+1,n+1);
 pcount = pcount(idxShuffle);
 while pcount(1) == 0                   %
   idxShuffle = randsample(n+1,n+1);    % somewhat embarassing: cannot fit 
   pcount = pcount(idxShuffle);         % if P(K=0) = 0. Need to fix this
 end                                    %
 
 output.spikes = sample_flat_model(pcount',size(output.spikes,2));
 
if ifApproxByBetaBinomial
    
  mu1 =  (0:n) * pcount;
  mu2 = (0:n).^2 * pcount ;
  Z = ( (n) * (mu2/mu1 - mu1 -1)) + mu1;
  a = (n * mu1 - mu2) / Z;
  b = (n - mu1) * (n - mu2/mu1) / Z; 
  lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';  
  logpcount = lognchoosek + betaln(a + (0:n), n + b - (0:n))' - betaln(a, b); 
  pcount = exp(logpcount);
  
  clear mu1 mu2 Z lognchoosek logpcount
end
  
if ifApproxByFlatIsing
     mu1 =  (0:n) * pcount;
     mu2 = (0:n).^2 * pcount ;
    [mu,rho,c]=meanvar_count_2_meancorr(mu1,mu2 - mu1^2,n);
    
    [hFI,JFI,pcount,~]=fit_flat_ising_model(mu,rho,n);
end  
pcountT = zeros(max(length(pcount)),length(Ns),length(Ts)); % P(K) at temperature T
varE = zeros(length(Ns),idxRep,length(Ts)); % storages for variance of energy/T^2 
cN = zeros(size(varE));             % and specific heat for each condition
varLogPr = zeros(length(Ts),length(Ns)); %
%varK = zeros(length(Ts),length(Ns));     % other stuff that might be interesting
h = zeros(length(Ts), length(Ns));       %
a = zeros(length(Ns),idxRep); 
b = zeros(length(Ns),idxRep); 

explogpcount = zeros(Ns(end)+1, length(Ns));
pcounts      = zeros(Ns(end)+1, length(Ns), idxRep);

lgnd = cell(length(Ns),1);   % legend for figure

lambdaT1 = -Inf * ones(n+1,length(Ns)); % lambda at T=1 for various N
lambdaT1(1,:) = 0;

% if ifRandomSubsampling
%   nSamples = 500000; % may take a few seconds... 
%   s = sample_flat_model(pcount',nSamples-nHack*n); % draw full samples from count distribution
%   s(:,end+1:end+n*nHack) = repmat(triu(ones(n)),1,nHack);
%   s = sparse(s);
%   clear output
% else

idxN = cell(length(Ns),1);
for j = 1:length(Ns)
    disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);
    if ifRandomSubsampling
        idxN{j} = zeros(Ns(j), idxRep);
        for t = 1:idxRep
            idxN{j}(:,t) = randsample(n,Ns(j));
        end
    else
        idxRep = 1;
        idxN{j} = idxSubsamples{Ns(j)/10}(:,t); % = (1:Ns(j))';
    end

    lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - gammaln(Ns(j)+1-(0:Ns(j))))';
        
    for t = 1:idxRep
        histCounts = full(sum(output.spikes(idxN{j}(:,t),:),1));
        hist_s = histc(histCounts,0:Ns(j))';
        hist_s = hist_s/sum(hist_s);
        
        if ifApproxByBetaBinomial
            
            mu1 =  (0:Ns(j)) * hist_s;
            mu2 = (0:Ns(j)).^2 * hist_s;
            Z = ( Ns(j) * (mu2/mu1 - mu1 -1)) + mu1;
            a(j,t) = (Ns(j) * mu1 - mu2) / Z;
            b(j,t) = (Ns(j) - mu1) * (Ns(j) - mu2/mu1) / Z;
            logpcount = lognchoosek + betaln(a(j,t) + (0:Ns(j)), Ns(j) + b(j,t) - (0:Ns(j)))' - betaln(a(j,t), b(j,t));
            explogpcount(1:Ns(j)+1, j) = exp(logpcount);            
            hist_s = exp(logpcount);
            
        end
        
        
        if ifApproxByFlatIsing
%            mu1 =  (0:Ns(j)) * hist_s;
%            mu2 = (0:Ns(j)).^2 * hist_s;
%            [mu,rho,~]=meanvar_count_2_meancorr(mu1,mu2 - mu1^2,Ns(j));            
%            [~,~,hist_s,~]=fit_flat_ising_model(mu,rho,Ns(j));
            [hist_s,~]=flat_ising_count_distrib(hFI,2*JFI,Ns(j));
        end
        
        pcounts(1:Ns(j)+1, j, t) = hist_s; % store
        
        
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
                    varE (j,t,i) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp.^2) - (pcountT(1:(Ns(j))+1,j,i)' * tmp)^2;
                    cN(j,t,i) = (1/Ns(j)) * varE(j,t,i); % specific heat, up to Boltzmann constant k_B
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

%clear E K T histCounts hist_s logpcount i j t idx lognchoosek logpcountT 


%% Preview figure for heat curves
%figure(42)
clrs = copper(18);

subplot(2,3,1:2),  
plot(0:length(pcount)-1, pcount, 'k', 'linewidth', 2.5)
title('shuffled P(K) of full population')
axis([0,length(pcount)-1, 0, 1.05*max(pcount)]);

subplot(2,3,4:5),
for i = 1:15, 
    plot(0:size(pcounts,1)-1, squeeze(pcounts(:,i,:)), 'color', clrs(19-i,:)); 
    hold on, 
end
title('shuffled P(K) for subpopulations')
axis([0,length(pcount)-1, 0, 1.05*max(pcounts(:))]);

subplot(2,3,[3,6]), 
for i = 1:15, 
    plot(Ts, squeeze(cN(i,:,:)), 'color', clrs(19-i,:)); 
    hold on, 
end
title('resulting heat curves')
axis([min(Ts), max(Ts), 0.95*min(cN(:)), 1.05*max(cN(:))]);


%%
clearvars -except Ts Ns cN pcounts pcount idxN idxShuffle

%save('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Figures/Fig5/figS_data_shuffle.mat', ...
%     'cN', 'pcount', 'pcounts', 'Ns', 'Ts', 'idxN', 'idxShuffle')