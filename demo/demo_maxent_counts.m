%fit maxent model with additional constraints on spike-count distribution

%% Run the code in which we need to construct the features by hand:
clear all
dim_x=10; %simulate 10 dimensional problem
Ns= 5000; %generate 1000 data-points;

h=randn(dim_x,1); %generate random bias terms;
J=randn(dim_x); J=triu(J,1)/sqrt(10); %generate random coupling terms:
%convert h and J to parameter lambda-- use convention that first n entries
%of lambda is a copy of h, and lower-triangular entries in J are ignored,
%and upper-triangular entries are first taken along the first row, then
%second row, etc:
lambda=hJ2lambda(h,J);

        
        
                
%set up second-order feature space for binary model:
[features,description,x]=setup_features_maxent(dim_x,2);
            
%calculate corresponding probabilities from ground-truth maximum entropy
%model:
[logPtrue,logZtrue,Ptrue, means_true]=logPMaxEnt(features,lambda);
logPtrueCheck=qIsing(x,h,J*2,exp(logZtrue),1);




[covo,meano]=wcov(all_states(dim_x),Ptrue);
 

%now, generate synthetic data from this distribution:
features_sampled=sample_discrete(features,Ptrue,max(Ns));

%calculate their means as we will need this as intput for the fitting
%procedure:
means_sampled=mean(features_sampled,1); clear features_sampled
                    

 
 %% We also have a function which does the bookkeeping internally, i.e. for whcih we never get to see the features: 

 [meano,covo]=meancov_2_features(means_sampled)
 means_check=meancov_2_features(meano,covo);


[h_checko,J_checko,logZ,logP, patterns]=fit_ising_model(meano,covo);

%%
figure
subplot(2,2,1)
plot(means_sampled,means_check,'.')
eqline
xlabel('Input features')
ylabel('features after converting twice')

subplot(2,2,2)
plot(h,h_checko,'.')
hold on
plot(J,J_checko,'g.')

eqline
xlabel('h J of first function')
ylabel('h J of second function')


pause(1)




