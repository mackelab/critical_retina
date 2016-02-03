% This script loads the spiking data of n simulated neurons, computes P(K) 
% from the spikes and then subsamples a new draw of patterns x (all neurons 
% are assumed identical w.r.t. firing rates and pairwise correlations) from
% P(K) to simulate activity at network sizes Ns(i)<=n, i = 1,...,length(Ns)

% As working with flat models is so simple, this file contains all the 
% flat_model related code lines.
 

clc
clear all
close all

% Load P(K) from simulation data: 
% The location on your hard-drive has to be changed accordingly!
load('../../data/Retina_model_data_new2.mat')
load('../../results/K_pairwise/idxSubsamples.mat')
disp('getting spatial subsampling order (may take a few seconds)')

n = size(cell_positions,2); % number of neurons in this simulation 

% ordering of cells (important for spatially structured subsampling
idxC = 1:size(cell_positions,2); % initialise as simple enumeration

% now get cell ordering along the length of the length of the retina patch
getSpatialSubsamplingOrder; % cleans up unused variables, 
                            % gives re-ordered indices idxC,
                            % also re-arranges output into struct variable,
                            % and starts plotting results

Ts = (0.8 : 0.0125 : 2)'; % range of temperatures to look at 

Ns = [20:20:300]; % range of network sizes to look at (max(Ns) should 
                  % be below size(pcount) = n+1, of course)

alphaOrT = 'T';     % decide whether to plainly vary the temperature T 
                    % or to do so while keeping the expected value of K
                    % constant by varying the 'alpha' in a related problem.
                    % Leave at 'T' if unsure, or consult Tkacik et al. 2014
                    % for details of varying said 'alpha'

ifApproxByBinomial     = false; % somewhat crude,                    
ifApproxByBetaBinomial = true;  % choose one
ifApproxByFlatIsing    = false; % (would not recommend choosing several)

% decide on how to generate subpopulations from full recording:
ifStructuredSubsampling = true; % replace random with 'spatially informed' 
% otherwise, we do uniform random subsampling and may draw many subpops
idxRep = 50; % number of subsets drawn for each n to estimate P(K) at n           
           
%--------------------------------------------------------------------------
% From here on, the proper algorithm starts (no need to change stuff)
%--------------------------------------------------------------------------

if  ifStructuredSubsampling 
 idxRep = 1; % does not make much sense to have this > 1 in this case
end

 pcount = sum(output.spikes,1)';
 pcount = histc(pcount,0:n);
 pcount = pcount/sum(pcount);
 
if ifApproxByBetaBinomial
  mu1 =  (0:n) * pcount;
  mu2 = (0:n).^2 * pcount ;
  Z = ( (n) * (mu2/mu1 - mu1 -1)) + mu1;
  a = (n * mu1 - mu2) / Z;
  b = (n - mu1) * (n - mu2/mu1) / Z; 
  lognchoosek = (gammaln(n+1) - gammaln((0:n)+1) - gammaln(n+1-(0:n)))';  
  logpcount = lognchoosek + betaln(a + (0:n), n+b-(0:n))'-betaln(a, b); 
  pcount = exp(logpcount);  
  clear mu1 mu2 Z lognchoosek logpcount  
end
  
if ifApproxByFlatIsing
     mu1 = (0:n) * pcount;
     mu2 = (0:n).^2 * pcount ;
    [mu,rho,c]=meanvar_count_2_meancorr(mu1,mu2 - mu1^2,n);    
    [hFI,JFI,pcount,~]=fit_flat_ising_model(mu,rho,n);
end  

% make space for computations
pcountT = zeros(max(length(pcount)),...
                length(Ns), ...
                length(Ts));  % P(K) at temperature T
varE = zeros(length(Ns),idxRep,length(Ts)); % storages for variance 
                                            % of energy/T^2 
varNcK  = zeros(length(Ns),idxRep,length(Ts)); 
varlPK  = zeros(length(Ns),idxRep,length(Ts)); 
covENcK = zeros(length(Ns),idxRep,length(Ts)); 
cN  = zeros(size(varE));             % and specific heat for each condition
cNx = zeros(size(varE));             
varLogPr = zeros(length(Ts),length(Ns)); %
h = zeros(length(Ts), length(Ns));       %
a = zeros(length(Ns),idxRep); 
b = zeros(length(Ns),idxRep); 

explogpcount = zeros(Ns(end)+1, length(Ns));
pcounts = explogpcount;

lgnd = cell(length(Ns),1);   % legend for figure

lambdaT1 = -Inf * ones(n+1,length(Ns)); % lambda at T=1 for various N
lambdaT1(1,:) = 0;

idxN = cell(length(Ns),1);
% for each population size ...
for j = 1:length(Ns)
    disp(['population size run #', num2str(j), ...
          ' out of ', num2str(length(Ns))]);
    if ifStructuredSubsampling
        idxN{j} = idxC(1:Ns(j));
    else
        idxN{j} = zeros(Ns(j), idxRep);
        for t = 1:idxRep
            idxN{j}(:,t) = randsample(n,Ns(j));
        end
    end

    lognchoosek = (gammaln(Ns(j)+1) - gammaln((0:Ns(j))+1) - ...
                   gammaln(Ns(j)+1-(0:Ns(j))))';
        
    % for each subsampled population ...
    for t = 1:idxRep
        histCounts = full(sum(output.spikes(idxN{j}(:,t),:),1));
        hist_s = histc(histCounts,0:Ns(j))';
        hist_s = hist_s/sum(hist_s); % data P(K)
        
        % check if to approximate P(K) with some analytic expression:
        if ifApproxByBinomial
            hist_s = binopdf(0:Ns(j), Ns(j), mean(histCounts)/Ns(j));
        end
        
        if ifApproxByBetaBinomial            
            mu1 =  (0:Ns(j)) * hist_s;
            mu2 = (0:Ns(j)).^2 * hist_s;
            Z = ( Ns(j) * (mu2/mu1 - mu1 -1)) + mu1;
            a(j,t) = (Ns(j) * mu1 - mu2) / Z;
            b(j,t) = (Ns(j) - mu1) * (Ns(j) - mu2/mu1) / Z;
            logpcount = lognchoosek + betaln(a(j,t) + ...
                        (0:Ns(j)), ...
                         Ns(j)+b(j,t)-(0:Ns(j)))'-betaln(a(j,t),b(j,t));
            explogpcount(1:Ns(j)+1, j) = exp(logpcount);            
            hist_s = exp(logpcount);
            
        end
                
        if ifApproxByFlatIsing
            [hist_s,~]=flat_ising_count_distrib(hFI,2*JFI,Ns(j));
        end
        
        pcounts(1:Ns(j)+1, j) = hist_s; % store
        
        
         % get lambda from count distr.
        lambdaT1((1:Ns(j))+1,j) = fit_flat_maxent_model(hist_s);
                
        for i = 1:length(Ts)
          T = Ts(i);
          lambdaT = lambdaT1(:,j)/T;
            
          switch alphaOrT
            case 'T'
              logpcountT = -Inf*ones(n+1,1);
              logpcountT(1:Ns(j)+1) = lambdaT(1:Ns(j)+1) + lognchoosek;
              pcountT(:,j,i) = exp(logpcountT); 
              % getting normaliser Z is easy here
              pcountT(:,j,i) = pcountT(:,j,i)/sum(pcountT(:,j,i)); 
              tmp = lambdaT(1:(Ns(j)+1));
              tmp(pcountT(1:(Ns(j))+1,j,i)==0) = 0;  % ensure 0 * Inf = 0
              varE(j,t,i) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp.^2) - ...
                             (pcountT(1:(Ns(j))+1,j,i)' * tmp)^2;
               % specific heat, up to Boltzmann constant k_B
              cN(j,t,i) = (1/Ns(j)) * varE(j,t,i);
              tmp1 = logpcountT(1:Ns(j)+1);
              tmp1(pcountT(1:Ns(j)+1,j,i)==0) = 0;   % ensure 0 * Inf = 0
              tmp2 = lognchoosek;
              varlPK(j,t,i) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp1.^2) - ...
                               (pcountT(1:(Ns(j))+1,j,i)' * tmp1)^2;
              varNcK(j,t,i) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp2.^2) - ...
                               (pcountT(1:(Ns(j))+1,j,i)' * tmp2)^2;
              covENcK(j,t,i)=  (pcountT(1:(Ns(j))+1,j,i)'*(tmp1.*tmp2)) ...
                              - (pcountT(1:(Ns(j))+1,j,i)' * tmp1) ...
                              * (pcountT(1:(Ns(j))+1,j,i)' * tmp2);
              cNx(j,t,i) = (varNcK(j,t,i) - 2 * covENcK(j,t,i) ...
                          + varlPK(j,t,i))/Ns(j);
            case 'alpha'
              [pcountT(1:(Ns(j)+1),j,i),~,models,h(i,j)]= ...
                                   flat_model_trace_temp(hist_s,1./T,true);
              varLogPr(i,j) = models.var_log_probs;
              pcountT(1:(Ns(j)+1),j,i) = pcountT(1:(Ns(j)+1),j,i) ...
                                           / sum(pcountT(1:(Ns(j)+1),j,i));
              % Compute variance of E via sampling from K-space:
              tmp = cumsum([0;squeeze(pcountT(1:(Ns(j)+1),j,i))]);
              [~,K] = histc(rand(nSamples,1), tmp);
              K = K-1;               % indexing...
              E = (lambdaT(K+1) + h(i,j) * K);
              varE(i,j,t) = var(E,t); % var(energy)/T^2
              % specific heat, up to Boltzmann constant k_B
              cN(i,j,t) = (1/Ns(j)) * varE(i,j,t); 
           end
        end
    end
    
end

clear E K T histCounts hist_s logpcount i j t idx lognchoosek logpcountT 

%% Preview figure for heat curves
figure(42)
subplot(1,6,(7:9)-6)
clrs = hsv(length(Ns));       % colors used by preview plot
ct = 1;
 for j = 1:length(Ns)
  for t = 1:idxRep
   plot(Ts, squeeze(cN(j,t,:)), 'color', clrs(j,:), 'marker', 'x', ...
        'linewidth', 1.5, 'markerSize', 6); hold on;
  end
%  lgnd{ct} = ['N = ', num2str(Ns(j))];
  ct = ct + 1;
 end
 axis([min(Ts)-0.1,max(Ts)+0.1,0.9*min(cN(:)),1.05*max(cN(:))])
 xlabel('temperature T'), ylabel('specific heat var[E]/(k_B T^2)')
 box off; set(gca, 'TickDir', 'out')
% legend(lgnd);
 hold off;
 title(['Heat curves'])

 clear pieakHi tmp 

subplot(1,6,(10:12)-6)
[~, idxT] = min(abs(Ts - 1));
plot(Ns, squeeze(cN(:,:,idxT))', 'bo-', 'linewidth', 2)
hold on
for t = 1:idxRep
 plot(Ns, squeeze(max(cN(:,t,:), [],3))', 'ko-', 'linewidth', 2)
end
legend('c(T=1)', 'max\{c(T)\}')
title('specific heat divergence vs. n')
xlabel('system size n')
ylabel('c(T)')
