%demo for running and testing 'flat' models, i.e. models which are
%completely specified by their spike-count distribution.

close all
clear all

%% Scenario 1: Sample from different flat DG models with same mean but different 
% correlation coefficients, and compare their rasters and
% count-distributions


N=100; %100 neurons
Nsamples=500; %1000 samples
mu=.1; %probability of spike of .1 per neuron
rhos=[0,.05,.1,.2,.4]; %pairwise correlations between neurons

for k=1:numel(rhos)
    [gamma(k),lambda(k),count_distrib{k},model(k)]=fit_flat_dg(mu,rhos(k),N);
    [s{k}]=sample_flat_model(count_distrib{k},Nsamples);
end


h(1)=figure

for k=1:5
subplot(2,3,k)
imagesc(-s{k}), colormap gray
xlabel('time')
ylabel('neurons');
title(['Correlation ', num2str(rhos(k))]);
end

subplot(2,3,6)
for k=1:5
   semilogy([0:N],model(k).count_distrib,'color',[0,0,k/5],'linewidth',2); 
   hold on
end
ylim([.0001,.2])





%% Scenario 2: Compare these models to Ising model


%% Scenario 3: Calculate entropy of these models as a function of population
% size, and compare to asymptotic results

%% 
