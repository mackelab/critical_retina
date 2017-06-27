% code for generation of summary figure on the criticality analysis (fig 0) 
% [figure 1 in current manuscript]

clear all
close all

load('../../results/K_pairwise_final/lambda_nat.mat')
load('../../results/K_pairwise_final/Efx_nat.mat')
ns = [2,6]; % pick population sizes (actual sizes are *10)
is = [7,1]; % pick (numbered) randomly drawn subpopulation 

%%

nSamples = [5000000, 500000]; 
burnIn = 1000; 
xSampled = cell(length(ns));
for j = 1:2
    
    n = 10 * ns(j);
    i = is(j);

    lambda = lambdaHatSim{n/10}(:,i);
    x = Efx{n/10}(:,i);

    EX = exp(lambda(1:n))./(1+exp(lambda(1:n))); 
    x0 = double(rand(n,1)<EX);

    tic;
    [xSampled{j},~,~] = maxEnt_gibbs_pair_C(nSamples(j), burnIn, lambda, x0, 'cluster');
    toc

    subplot(1,2,j)
    plot(x, xSampled{j}, 'x')
    line([min(x), max(x)], [min(x), max(x)], 'color', 'k')
end
%save('fig1_data.mat', 'ns', 'is', 'xSampled')
%%
clc
clear all
load('../../results/K_pairwise_final/lambda_nat.mat')
load('../../results/K_pairwise_final/Efx_nat.mat')
load('fig1_data.mat')

clrs = [254,153,41;236,112,20;204,76,2;153,52,4;102,37,6;0,0,0]/255;

figure; 
for j = 2:-1:1,
    
    n = 10 * ns(j);
    i = is(j);    
    
    fDescrJ = nchoosek(1:n,2)'; 
    covsx = Efx{n/10}((n+1):(n*(n+1)/2),i) - ... % covariance as computed
           (Efx{n/10}(fDescrJ(1, :),i).* ...
            Efx{n/10}(fDescrJ(2, :),i)); % from Efx    for j = js, 
    varsx = Efx{n/10}(1:n,i) .* (1 - Efx{n/10}(1:n,i));
    corsx = covsx ./(sqrt(varsx(fDescrJ(1,:))).*sqrt(varsx(fDescrJ(2,:))));
    Efy = xSampled{j};
    covsy = Efy((n+1):(n*(n+1)/2)) - ... % covariance as computed
           (Efy(fDescrJ(1, :)).* Efy(fDescrJ(2, :))); % from Efx    
    varsy = Efy(1:n) .* (1 - Efy(1:n));
    corsy = covsy ./(sqrt(varsy(fDescrJ(1,:))).*sqrt(varsy(fDescrJ(2,:))));
    subplot(3,1,1)
    if j == 2
        line([0, 5], [0, 5], 'color', 'k')
        hold on
    end
    plot(50 * Efx{n/10}(1:n,i), 50 * Efy(1:n), 'o', 'color', clrs(n/20,:), ...
        'markersize', 3, 'markerFaceColor', clrs(n/20,:))
    axis square
    hold on;
    axis([0.015, 0.08, 0.015, 0.08] * 50)
    set(gca, 'XTick', [1:4])
    set(gca, 'YTick', [1:4])
    set(gca, 'TickDir', 'out')
    xlabel('Firing rate [Hz]')
    ylabel('Firing rate [Hz]')
    box off
    
    subplot(3,1,2)
    if j == 2
        line([-0.1, 0.6], [-0.1, 0.6], 'color', 'k')
        hold on
    end
    plot(corsx, corsy, 'o', 'color', clrs(n/20,:), ..., 
        'markersize', 3, 'markerFaceColor', clrs(n/20,:))
    axis square
    axis([-0.05, 0.6, -0.05, 0.6])
    hold on;
    set(gca, 'XTick', [0, 0.2, 0.4])
    set(gca, 'YTick', [0, 0.2, 0.4])
    set(gca, 'TickDir', 'out')
    xlabel('correlation')
    ylabel('correlation')
    box off

    subplot(3,1,3)        
    semilogy(0:n, Efx{n/10}(end-n:end,i), '.-', 'color', 0.7 * [1,1,1])
    hold on
    semilogy(0:n, Efy(end-n:end), '.-', 'color', clrs(n/20,:));
    axis square
    axis([0, 33, 0, 1])
    axis autoy
    hold on;
    xlabel('population spike-count')
    ylabel('probability')
    set(gca, 'XTick', [0, 10, 20, 30])
    set(gca, 'YTick', [10^-5, 10^-3, 10^-1])
    set(gca, 'TickDir', 'out')
    box off
    
end
