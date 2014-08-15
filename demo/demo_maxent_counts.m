%fit maxent model with additional constraints on spike-count distribution

%% First generate some samples from an Ising model: 
clear all
close all

dim_x=8; %simulate 10 dimensional problem
Ns= 100000; %generate 1000 data-points;

h=randn(dim_x,1)-2; %generate random bias terms;
J=randn(dim_x); J=triu(J,1)/sqrt(10); 
lambda=hJ2lambda(h,J);

        
%set up second-order feature space for binary model:
[features,description,x]=setup_features_maxent(dim_x,2);
            
%calculate corresponding probabilities from ground-truth maximum entropy
%model:
[logPtrue,logZtrue,Ptrue, means_true]=logPMaxEnt(features,lambda);
%to check, also 
logPtrueCheck=qIsing(x,h,J,exp(logZtrue),1);

[cov_true,mean_true]=wcov(x,Ptrue);
 
%now, generate synthetic data from this distribution:
x_sampled=sample_discrete(x,Ptrue,max(Ns));
%calculate  feature-means for ising_count model (not that this is a
%horribly memory-inefficient way of doing it, but I dont care at the moment)

[features_sampled]=setup_features_maxent(x_sampled,'ising_count');
[features_count]=setup_features_maxent(x,'ising_count');
means_sampled=mean(features_sampled,1); clear features_sampled
 
count_histogram_true=ps_2_count_distrib(Ptrue,x);
count_histogram_sampled=ps_2_count_distrib(ones(Ns,1)/Ns,x_sampled);

fitoptions.optTol=1e-100; 
fitoptions.progTol=1e-100; 
fitoptions.display='off';
fitoptions.MaxIter=3000;

[lambda_learned,logZlearned, logPlearned, means_learned,output]=fit_maxent_linear(features_count,means_sampled,fitoptions);
[lambda_learned_ising,logZlearned_ising, logPlearned_ising, means_learned_ising,output_ising]=fit_maxent_linear(features_count(:,1:end-dim_x),means_sampled(1:end-dim_x),fitoptions);
 
[logPlearned]=logPMaxEnt(features_count,lambda_learned);
[logPlearned_ising]=logPMaxEnt(features,lambda_learned_ising);
count_histogram_learned=ps_2_count_distrib(exp(logPlearned),x);
count_histogram_ising=ps_2_count_distrib(exp(logPlearned_ising),x);



%%
figure
subplot(2,2,1)
plot(0:dim_x,count_histogram_true,'-')
hold on
plot(0:dim_x,count_histogram_sampled,'r-')
plot(1:dim_x,means_sampled(end-dim_x+1:end),'o')

xlabel('Counts')
ylabel('Frequency')

subplot(2,2,2)
plot(logPtrue,logPlearned,'b.')
hold on
plot(logPtrue,logPlearned_ising,'r.')
eqline


subplot(2,2,3)
plot(lambda,'.')
hold on
plot(lambda_learned,'g.')
plot(lambda_learned_ising,'k.')

subplot(2,2,4)
semilogy(0:dim_x,count_histogram_true,'-')
hold on
plot(0:dim_x,count_histogram_learned,'r-')
plot(0:dim_x,count_histogram_ising,'o')

