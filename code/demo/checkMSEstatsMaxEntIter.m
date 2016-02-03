
%% 

 load('C:\Users\Loki\Desktop\Criticality\Code_critical_retina\Heat curves\EfxSimData.mat')
 load('C:\Users\Loki\Desktop\Criticality\Code_critical_retina\Heat curves\EfxSimData.mat')
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
eps = [0.01, 0.01, 0.05]; 
 clrs = hsv(3);
 load('C:\Users\Loki\Desktop\Criticality\Code_critical_retina\Heat curves\EfxSimData.mat')
 fig1 = figure('units','normalized','position',[0,0,0.6,0.9]);
 hndl = [];
for i = 1:10
 delete(hndl)
n = 10*i;
load(['C:\Users\Loki\Desktop\tmp\AssemblyLine\hJVn', num2str(n), '.mat'])

fDescrJ = nchoosek(1:n,2)'; 
for j = 1:length(fD)
 covs.x = Efx{i}((n+1):(n*(n+1)/2),1) - ...              % covariance as computed
        (Efx{i}(fDescrJ(1, :),1).* Efx{i}(fDescrJ(2, :),1)); % from Efx    
    
 subplot(231), plot(0:size(fD{j}.MSEperc,2)-2, fD{j}.MSEperc(1,2:end), 'color', clrs(j,:));   title('E[X] for (h,J,V) model')
 ylabel('MSE^{0.5} / (1/n \lambda^2)^{0.5}')
 box off           
 hold on
 axis([0, size(fD{j}.MSEperc,2)-1, 0, 0.2])
 line([0,size(fD{j}.MSEperc,2)], [eps(2), eps(2)], 'color', 'k')
 subplot(232), plot(0:size(fD{j}.MSEperc,2)-2, fD{j}.MSEperc(3,2:end), 'color', clrs(j,:));   title('P(K) for (h,J,V) model')
 box off           
 hold on           
 axis([0, size(fD{j}.MSEperc,2)-1, 0, 0.2])
 line([0,size(fD{j}.MSEperc,2)], [eps(3), eps(3)], 'color', 'k')
 subplot(233), plot(0:size(fD{j}.MSEperc,2)-2, fD{j}.MSEperc(4,2:end), 'color', clrs(j,:)); title('cov(X) for (h,J,V) model')
 box off           
 hold on
 axis([0, size(fD{j}.MSEperc,2), 0, 0.5])
 line([0,size(fD{j}.MSEperc,2)], [eps(2), eps(2)], 'color', 'k')

end


 load(['C:\Users\Loki\Desktop\tmp\AssemblyLine\hVn', num2str(n), '.mat'])
for j = 1:length(fD) 
 subplot(234), plot(0:size(fD{j}.MSEperc,2)-2, fD{j}.MSEperc(1,2:end), 'color', clrs(j,:));    title('E[X] for (h,V) model')
 xlabel('#iter')
 ylabel('MSE^{0.5} / (1/n \lambda^2)^{0.5}')
 box off           
 hold on
 axis([0, size(fD{j}.MSEperc,2)-1, 0, 0.2])
 line([0,size(fD{j}.MSEperc,2)], [eps(1), eps(1)], 'color', 'k')
 subplot(235), plot(0:size(fD{j}.MSEperc,2)-2, fD{j}.MSEperc(3,2:end), 'color', clrs(j,:));    title('P(K) for (h,V) model')
 xlabel('#iter')
 box off           
 hold on
 axis([0, size(fD{j}.MSEperc,2)-1, 0, 0.2])
 line([0,size(fD{j}.MSEperc,2)], [eps(3), eps(3)], 'color', 'k')
 subplot(236), plot(0:size(fD{j}.MSEperc,2)-2, fD{j}.MSEperc(4,2:end), 'color', clrs(j,:));  title('cov(X) for (h,V) model')
 xlabel('#iter')
 box off           
 hold on
 axis([0, size(fD{j}.MSEperc,2)-1, 0, 0.5])
 line([0,size(fD{j}.MSEperc,2)], [eps(2), eps(2)], 'color', 'k')
end
 
hndl = annotation('textbox', [0 0.9 1 0.1], ...
      'String', ['Error metrics for n = ', num2str(n), ', horizontal lines give chosen threshold'], ...
      'EdgeColor', 'none', ...
      'HorizontalAlignment', 'center');

 for k = 1:6
  subplot(2,3,k), hold off;
  set(gca, 'XTick', get(gca, 'XTick'))
  tcklbl = get(gca, 'XTickLabel');
  set(gca, 'XTickLabel', 25*str2num(tcklbl)) 
 end
 
 export_fig(['C:\Users\Loki\Desktop\tmp\resAnal',num2str(n),'.pdf'])


 pause;

end 