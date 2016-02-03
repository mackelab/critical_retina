% This script loads the spiking data of n simulated neurons, computes P(K) 
% from the spikes and then subsamples a new draw of patterns x (all neurons 
% are assumed identical w.r.t. firing rates and pairwise correlations) from
% P(K) to simulate activity at network sizes Ns(i)<=n, i = 1,...,length(Ns)

clc
clearvars -except Output Res
%close all

% Load P(K) from simulation data: 
% The location on your hard-drive has to be changed accordingly!

% n = 316 cells (small simulation), natural stimuli 
%load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_data_new2.mat')

% n = 316 cells (small simulation), checkerboard stimuli
%load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_data_checkerboard_alt2.mat')

% n = 316 cells (small simulation), full-field flicker stimuli
load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_data_flicker_alt.mat')

% n = 6000 cells (big simulation), natural stimuli 
%load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_gaussian_neu_offset_029.mat')

% n = 6000 cells (big simulation), checkerboard stimuli 
%load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_gaussian_checker_offset_029.mat')

% n = 6000 cells (big simulation), full-field flicker stimuli 
%load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Model Retina/Retina_model_gaussian_FF_offset_029.mat')

disp('getting spatial subsampling order (may take a few seconds)')

%idxC = abs(cell_positions(2,:)) < 50;
%idxC = abs(cell_positions(1,:)) < 400;
idxC = 1:size(cell_positions,2);
output = output(:, :, idxC);
cell_positions = cell_positions(:,idxC);

getSpatialSubsamplingOrder; % cleans up unused variables, 
                            % gives re-ordered indices idxC,
                            % also re-arranges output into struct variable,
                            % and starts plotting results
%output.spkCorrs = [];

clearvars -except idxC output cell_positions m Output Res
%%
n = size(cell_positions,2); % number of neurons in this simulation 
if n > 316 % big simulation
    n = 6000; % limit to first 6000 cells from the left
end

Ts = (0.9 : 0.00125 : 1.4)'; % range of temperatures to look at 

%Ns = [15,30,60,120,240];% range of network sizes to look at (max(Ns) should 
Ns = [20:4:n];          % be below size(pcount) = n+1, of course)
                        
alphaOrT = 'T';     % decide whether to plainly vary the temperature T 
                    % or to do so while keeping the expected value of K
                    % constant by varying the 'alpha' in a related problem.
                    % Leave at 'T' if unsure, or consult Tkacik et al. 2014
                    % for details of varying said 'alpha'

ifApproxByBetaBinomial = true; 

ifRandomSubsampling    = true; % replace 'spatially informed' with random
idxRep = 50; % number of subsets drawn for each n to estimate P(K) at n

    

if ifRandomSubsampling % overwrite cell re-ordering with default 1:n
 idxC =  1:length(idxC);                           % all cells
end

% reduce considered area of retina patch and/or apply re-ordering of cells
cell_positions = cell_positions(:, idxC); output.spikes = output.spikes(idxC,:);
%idxC = 1:length(idxC); % in case of repeated evaluations of this part, don't want to 
%                       % re-shuffle the cells over and over again.

nHack = 1; % number of guaranteed samples for each K>0 to sneak into 
           % the data (makes sure P(K)>0)

%--------------------------------------------------------------------------
if  ~ifRandomSubsampling 
 idxRep = 1; % does not make much sense to have this > 1 in this case
end

% From here on, the proper algorithm starts (no need to change stuff)
  
pcountT = zeros(n+1,length(Ns),length(Ts)); % P(K) at temperature T
varE = zeros(length(Ns),idxRep,length(Ts)); % storages for variance of energy/T^2 
cN = zeros(size(varE));             % and specific heat for each condition
varLogPr = zeros(length(Ts),length(Ns)); %
h = zeros(length(Ts), length(Ns));       %
a = zeros(length(Ns), idxRep); 
b = zeros(length(Ns), idxRep); 

explogpcount = zeros(Ns(end)+1, length(Ns));
pcounts = explogpcount;

lgnd = cell(length(Ns),1);   % legend for figure

lambdaT1 = -Inf * ones(n+1,length(Ns)); % lambda at T=1 for various N
lambdaT1(1,:) = 0;


idxN = cell(length(Ns),1);

for t = 1:idxRep
  idxN{end}(:,t) = randsample(n,Ns(end));
end

for j = 1:length(Ns)
    disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);
    if ifRandomSubsampling
        if j < length(Ns)
         idxN{j} = zeros(Ns(j), idxRep);
        end
        for t = 1:idxRep
            idxN{j}(:,t) = idxN{end}(1:Ns(j),t);
        end
        
    else
        idxRep = 1;
        idxN{j} = (1:Ns(j))';
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
            if Z < 0 || a(j,t) < 0 || b(j,t) < 0
              [s,~]=sample_flat_model(hist_s', 10000);
              [a(j,t),b(j,t)] = fitBetabinML(sum(s,1), Ns(j), 10000);
              disp(['method of moments failed for n = ', num2str(Ns(j))])
            end
            logpcount = lognchoosek + betaln(a(j,t) + (0:Ns(j)), Ns(j) + b(j,t) - (0:Ns(j)))' - betaln(a(j,t), b(j,t));
            explogpcount(1:Ns(j)+1, j) = exp(logpcount);            
            hist_s = exp(logpcount);
        end
        
        
        pcounts(1:Ns(j)+1, j) = hist_s; % store
        
        
        lambdaT1((1:Ns(j))+1,j) = fit_flat_maxent_model(hist_s); % get lambda from count distr.                
        
        for i = 1:length(Ts)
            T = Ts(i);
            lambdaT = lambdaT1(:,j)/T;
            
                    logpcountT = -Inf*ones(n+1,1);
                    logpcountT(1:Ns(j)+1) = lambdaT(1:Ns(j)+1) + lognchoosek;
                    pcountT(:,j,i) = exp(logpcountT);
                    pcountT(:,j,i) = pcountT(:,j,i)/sum(pcountT(:,j,i)); % getting Z is easy here
                    h(i,j) = 0;
                    tmp = lambdaT(1:(Ns(j)+1));
                    tmp(pcountT(1:(Ns(j))+1,j,i)==0) = 0;
                    varE (j,t,i) =  (pcountT(1:(Ns(j))+1,j,i)' * tmp.^2) - (pcountT(1:(Ns(j))+1,j,i)' * tmp)^2;
                    cN(j,t,i) = (1/Ns(j)) * varE(j,t,i); % specific heat, up to Boltzmann constant k_B
        end
    end
    
end

clear E K T histCounts hist_s logpcount i j t idx lognchoosek logpcountT 


%% Preview figure for heat curves
figure(42)
[~, idxT] = min(abs(Ts - 1));

clrs = [254,153,41;236,112,20;204,76,2;153,52,4;102,37,6;0,0,0]/255;

plot(-1,-1, 'color', clrs(2,:), 'linewidth', 2.5)
hold on
plot(-1,-1, 'color', clrs(4,:), 'linewidth', 2.5)

tmp = zeros(length(Ns), idxRep);
for t = 1:idxRep
 tmp(:,t) = squeeze(max(cN(:,t,:), [],3));
end
if idxRep > 1
    h = area(Ns, ...
         [ mean(tmp,2)-  std(tmp,0,2), ...
                       2*std(tmp,0,2)]);
    h(1).FaceColor = [1,1,1];
    h(1).EdgeColor = 'none';
    h(2).FaceColor = [clrs(4,:)];
    h(2).EdgeColor = 'none';
end
plot(Ns, mean(tmp,2), '-',  'color', 'k', 'linewidth', 2)
if idxRep < 2
    plot(Ns, mean(tmp,2), 's',  'color', clrs(4,:), 'linewidth', 2)
end
if idxRep > 1
    h = area(Ns, ...
         [ mean(squeeze(cN(:,:,idxT)),2)-  std(squeeze(cN(:,:,idxT)),0,2), ...
                                         2*std(squeeze(cN(:,:,idxT)),0,2)]);
    h(1).FaceColor = [1,1,1];
    h(1).EdgeColor = 'none';
    h(2).FaceColor = [clrs(2,:)];
    h(2).EdgeColor = 'none';
end
plot(Ns, mean(squeeze(cN(:,:,idxT)),2), '-', 'color', 'k', 'linewidth', 2)
if idxRep < 2
    plot(Ns, mean(squeeze(cN(:,:,idxT)),2), 's',  'color', clrs(2,:), 'linewidth', 2)
end
               
%line([Ns(1), Ns(end)], [cN(idxT, 1), cN(idxT, end)], 'color', 'b')     % reference lines (will be very similar
%line([Ns(1), Ns(end)], [max(cN(:, 1)), max(cN(:, end))], 'color', 'k') % to the results from random subsampling)
legend('c(T=1)', 'max\{c(T)\}')
title('specific heat divergence against system size n (should saturate)')
xlabel('system size n')
ylabel('c(T)')

axis([min(Ns)-1, max(Ns)+1, 0, 1]); axis autoy
legend('c(T=1)', 'max\{c(T)\}', 'Location', 'Northwest')
legend boxoff
set(gca, 'TickDir', 'out')
box off

