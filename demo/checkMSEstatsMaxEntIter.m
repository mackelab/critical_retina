
%% 

 load('C:\Users\Loki\Desktop\tmp\MSEsFirstSimRuns')
 load('EfxSimData.mat')
 figure;

for i = 3:10
n = 10*i;    
fDescrJ = nchoosek(1:n,2)'; 
covs.x = Efx{i}((n+1):(n*(n+1)/2),1) - ...              % covariance as computed
        (Efx{i}(fDescrJ(1, :),1).* Efx{i}(fDescrJ(2, :),1)); % from Efx    
    
 subplot(231), plot(sqrt(MSEhJV{i}.h)/sqrt(mean(Efx{i}(1:n,1).^2)));   title('E[X] for (h,J,V) model')
 ylabel('log(MSE')
 box off           
 axis([0, length(MSEhJV{i}.h), 0, 0.1])
 subplot(232), plot(sqrt(MSEhJV{i}.V)/sqrt(mean(Efx{i}(end-n:end,1).^2)));   title('P(K) for (h,J,V) model')
 box off           
 axis([0, length(MSEhJV{i}.V), 0, 0.1])
 subplot(233), plot(sqrt(MSEhJV{i}.cov)/sqrt(mean(abs(covs.x).^2))); title('cov(X) for (h,J,V) model')
 box off
 axis([0, length(MSEhJV{i}.cov), 0, 0.1])

 subplot(234), plot(sqrt(MSEhV{i}.h)/sqrt(mean(Efx{i}(1:n,1).^2)));    title('E[X] for (h,V) model')
 xlabel('#iter')
 ylabel('log(MSE')
 box off
 axis([0, length(MSEhV{i}.h), 0, 0.1])
 subplot(235), plot(sqrt(MSEhV{i}.V)/sqrt(mean(Efx{i}(end-n:end,1).^2)));    title('P(K) for (h,V) model')
 xlabel('#iter')
 box off
 axis([0, length(MSEhV{i}.V), 0, 0.1])
 subplot(236), plot(sqrt(MSEhJV{i}.cov)/sqrt(mean((covs.x).^2)));  title('cov(X) for (h,V) model')
 xlabel('#iter')
 box off
 axis([0, length(MSEhV{i}.cov), 0, 0.1])
 
 pause
end

%%


fDescrJ = nchoosek(1:n,2)'; 
covs.x = Efx((n+1):(n*(n+1)/2),1) - ...              % covariance as computed
        (Efx(fDescrJ(1, :),1).* Efx(fDescrJ(2, :),1)); % from Efx    
    
 subplot(231), plot(sqrt(MSEhJV.h)/sqrt(mean(Efx(1:n,1).^2)));   title('E[X] for (h,J,V) model')
 ylabel('log(MSE')
 box off           
 axis([0, length(MSEhJV.h), 0, 0.1])
 subplot(232), plot(sqrt(MSEhJV.V)/sqrt(mean(Efx(end-n:end,1).^2)));   title('P(K) for (h,J,V) model')
 box off           
 axis([0, length(MSEhJV.V), 0, 0.1])
 subplot(233), plot(sqrt(MSEhJV.cov)/sqrt(mean(abs(covs.x).^2))); title('cov(X) for (h,J,V) model')
 box off
 axis([0, length(MSEhJV.cov), 0, 0.1])

 subplot(234), plot(sqrt(MSEhV.h)/sqrt(mean(Efx(1:n,1).^2)));    title('E[X] for (h,V) model')
 xlabel('#iter')
 ylabel('log(MSE')
 box off
 axis([0, length(MSEhV.h), 0, 0.1])
 subplot(235), plot(sqrt(MSEhV.V)/sqrt(mean(Efx(end-n:end,1).^2)));    title('P(K) for (h,V) model')
 xlabel('#iter')
 box off
 axis([0, length(MSEhV.V), 0, 0.1])
 subplot(236), plot(sqrt(MSEhJV.cov)/sqrt(mean((covs.x).^2)));  title('cov(X) for (h,V) model')
 xlabel('#iter')
 box off
 axis([0, length(MSEhV.cov), 0, 0.1])
 