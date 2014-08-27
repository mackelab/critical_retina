%clear all
%close all

d=10; %simulate 10 dimensional problem
n= 10000; %generate 1000 data-points;
model = 'ising_count_l_0';
newLambda = true;

if newLambda
 h=randn(d,1)-1; %generate random bias terms;
 J=randn(d); J=triu(J,1)/sqrt(d); 
 lambda=hJ2lambda(h,J);
 switch model
     case 'ising_count_l_0'
         L=randn(d+1,1)/sqrt(d);
     case 'ising_count'
         L=randn(d,1)/sqrt(d);
     case 'ising'
         L = [];
 end
 lambda = [lambda;L];
end

[features,description,x]=setup_features_maxent(d,model);
[logPtrue,logZtrue,Ptrue, means_true]=logPMaxEnt(features,lambda);

[cov_true,mean_true]=wcov(x,Ptrue);
%now, generate synthetic data from this distribution:§
x_sampled=sample_discrete(x,Ptrue,max(n));
fx_sampled = setup_features_maxent(x_sampled, model);

figure;  
subplot(1,3,1), plot(mean(x_sampled,1), 'ko-', 'linewidth', 3)
axis([0.5, d+0.5, 0, 1]);
box off
title('Getting the marginals:  P(x_i), i = 1, ..., d')
set(gca, 'TickDir', 'out')

subplot(1,3,2), plot(mean(fx_sampled,1), 'ko-', 'linewidth', 3)
axis([0.5, size(fx_sampled,2)+0.5, 0, 1]);
box off
title('Getting the moments:  E[f_i(x)], i = 1, ..., d')
set(gca, 'TickDir', 'out')
h = histc(sum(x_sampled,2),0:d)/n;
subplot(1,3,3), plot(0:d, h, 'ko-', 'linewidth', 3)
box off
title('Getting the counts:  P[K], K = 0, ..., d')
set(gca, 'TickDir', 'out')
axis([-0.5, d+0.5, 0, 1.1*max(h)]);


nRestarts = 4;
clrs = hsv(nRestarts);
for i = 1:nRestarts
disp(['Restart number ', num2str(i), ' of ', num2str(nRestarts)])
nSamples = 100000; burnIn = 50000; x0 = x_sampled(end,:)'; 
[xSampled,p1s,p0s,ks] = maxEnt_gibbs(nSamples, burnIn, -lambda, x0, model);
fxSampled = setup_features_maxent(xSampled', model);
figure(1); 
subplot(1,3,1), hold on; plot(mean( xSampled,2), 'o-', 'color', clrs(i,:), 'linewidth', 2)
subplot(1,3,2), hold on; plot(mean(fxSampled,1), 'o-', 'color', clrs(i,:), 'linewidth', 2)
subplot(1,3,3), hold on; plot(0:d, histc(sum(xSampled,1),0:d)/nSamples, 'o-', 'color', clrs(i,:), 'linewidth', 2)
%figure(2); 
% for k = 1:d
%     id = ks==k;
%     subplot(d,2,2*k-1), hold on; plot(histc(p1s(id),0:0.1:1.5), 'color', clrs(i,:)); 
%     subplot(d,2,2*k), hold on; plot(0:0.05:1, histc(p1s(id)./(p0s(id)+p1s(id)), 0:0.05:1), 'color', clrs(i,:))
% end
end

% Sanity check: compare count distribution to one that is expected if one
% actually had fitted an independent model (J_ij = 0 for all i,j)
means = mean( xSampled,2);
tmp = zeros(nSamples,1);
for t = 1:nSamples % for each samples ...
 for k = 1:d % ... draw d times independently from E[X_k = 1] ...
  tmp(t) = tmp(t) + (rand(1)<means(k)); % ... and add all up
 end
end
plot(0:d, histc(tmp, 0:d)/nSamples, 'kx--');

switch model
    case 'ising'
        model_name = 'Ising';
    case 'ising_count'
        model_name = 'Ising + V(K)';
    case 'ising_count_l_0'
        model_name = 'Ising + V(K)';
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Gibbs sampling results for ', model_name, ', d = ', ...
                num2str(d), ', #Samples = ', num2str(nSamples), ...
                ', #burnIn = ', num2str(burnIn)], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
