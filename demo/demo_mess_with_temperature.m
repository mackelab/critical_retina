%
close all
clear all

N=200; %100 neurons
Nsamples=100; %1000 samples
mu=.3; %probability of spike per neuron
rho=.1; %pairwise correlation between neurons

betas=10.^[-.5:.01:.5];
set(gcf,'DefaultAxesColorOrder',jet(round(numel(betas(1:10:end)))));

[gamma,lambda,DG_probs,dg_model]=fit_flat_dg(mu,rho,N);
[s]=sample_flat_model(DG_probs,Nsamples);
    
[lambda,maxent_probs,maxent_model]=fit_flat_maxent_model(dg_model.count_distrib);
[maxent_s]=sample_flat_model(maxent_probs,Nsamples);

[count_distribs,Zs,temp_models]=flat_model_trace_temp(dg_model.count_distrib,betas);





subplot(2,4,1)
semilogy([0:N],dg_model.count_distrib,'linewidth',2); 
title('Count distribution, original model')   

subplot(2,4,2)
semilogx(betas,[temp_models.mean])
title('Mean as function of beta')   

subplot(2,4,3)
semilogx(betas,[temp_models.corr])
title('Corr as function of beta') 

subplot(2,4,4)
semilogx(betas,[temp_models.entropy])
title('Entropy as function of beta') 

subplot(2,4,5)
semilogx(betas,[temp_models.var_log_probs])
title('Var(log probs) as function of beta') 

subplot(2,4,6)
semilogy([0:N],vertcat(temp_models(1:10:end).count_distrib))
hold on
semilogy([0:N],dg_model.count_distrib,'linewidth',3); 
ylim([1e-30,1])