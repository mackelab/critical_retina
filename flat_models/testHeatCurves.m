Ts = (0.5 : 0.05 : 4)';
Ns = 10:10:100;
%load('C:\Users\Loki\Desktop\Criticality\Code\Package\Data\SalamanderDataStatistics.mat');
%pcount = data.populationAcitivityCounts_valuesForKfrom0to40'; % get true count distribution
load('C:\Users\Loki\Desktop\Criticality\Data\tkacik2014N100logPK.mat');
pcount = zeros(101,1);
pcount(1:46) = 10.^(tkacik2014N100logPK(:,2)');
pcountT = zeros(max(length(pcount)),length(Ns));
varE = zeros(length(Ts),length(Ns));
CN = zeros(size(varE));

n = length(pcount)-1; % number of neurons

figure(1);
lgnd = cell(length(Ns),1);
textdx = 0.05; textdy = 0.1;
% getting lambda the awkward way (via full distribution over n-dimensional neural states...)
nSamples = 10000;
[s]=sample_flat_model(pcount',nSamples); % draw full samples from count distribution
for j = 1:length(Ns)
  disp(['population size run #', num2str(j), ' out of ', num2str(length(Ns))]);
  idxN = randsample(n,Ns(j));
  hist_s = zeros(Ns(j)+1,1); 
  histCounts = sum(s(idxN,:),1);
  for l = 0:Ns(j);
    hist_s(l+1) = mean(histCounts==l);
  end
  lambdaT1 = [0;fit_flat_maxent_model(hist_s)]; % get lambda from count distr.

 for i = 1:length(Ts)
    disp(['- temperature run #', num2str(i), ' out of ', num2str(length(Ts))]);
    T = Ts(i);
    lambdaT = lambdaT1/T;
    for k = 0:Ns(j)
        pcountT(k+1,j) =  nchoosek(Ns(j),k) * exp(lambdaT(k+1));
    end
    pcountT(:,j) = pcountT(:,j) / sum(pcountT(:,j));
    si=sample_flat_model(pcountT(1:(Ns(j)+1),j)',nSamples);
    E=lambdaT(sum(si,1)+1);
    varE(i,j) = var(E);
 end
 CN(:,j) = varE(:,j)./(T.^2);
end

for j = 1:length(Ns)
 plot(Ts, CN(:,j), 'color', 'k', 'marker', 'o', 'linewidth', 1.5, 'markerSize', 1.5); hold on;
end
boxl = 0.75*max(Ts); boxr = max(Ts); boxu = 2/3*max(CN(:)); boxd = 1/6*max(CN(:));
for j = 1:length(Ns)
 [peakHi,peakLoc] = max(CN(:,j));
 line([Ts(peakLoc),Ts(peakLoc)],[peakHi,peakHi+0.25],'color','g','linewidth',1.5);
 line([Ts(peakLoc),Ts(peakLoc)+0.02],[peakHi,peakHi+0.1],'color','g','linewidth',1.5);
 line([Ts(peakLoc),Ts(peakLoc)-0.02],[peakHi,peakHi+0.1],'color','g','linewidth',1.5);
 text(Ts(peakLoc)+textdx, peakHi+textdy, ...
      ['T = ', num2str(Ts(peakLoc)), ', C/N = ', num2str(peakHi), ', C/N^2 = ', num2str(peakHi/Ns(j))], ...
       'color', 'g');
 lgnd{j} = ['N = ', num2str(Ns(j))];
 
 plot(boxl+0.1*(boxr-boxl)+0.8*(Ns(j))/(max(Ns))*(boxr-boxl), boxd+0.1*(boxu-boxd)+0.8*peakHi/max(CN(:))*(boxu-boxd), 'g*');
end
text(boxl+0.49*(boxr-boxl), boxd-0.1, 'N');
text(boxl+0.17*(boxr-boxl), boxu+0.1, 'Dependence of peak C/N(T) on N');
text(boxl-0.1, boxd+0.5*(boxu-boxd), 'C/N');
line([1,1],[0,1.05*max(varE(:))], 'color', 'r', 'linewidth', 1.5)
line([boxl,boxl],[boxu,boxd], 'color', 'k');
line([boxr,boxr],[boxu,boxd], 'color', 'k');
line([boxl,boxr],[boxu,boxu], 'color', 'k');
line([boxl,boxr],[boxd,boxd], 'color', 'k');
line([boxl,boxr],[boxd,boxu], 'color', 'b');
axis([min(Ts)-0.1,max(Ts)+0.1,0.9*min(CN(:)),1.05*max(CN(:)+textdy)])
xlabel('temperature T'), ylabel('specific heat var[E]/(k_B T^2)')
box off; set(gca, 'TickDir', 'out')
legend(lgnd);
hold off;
