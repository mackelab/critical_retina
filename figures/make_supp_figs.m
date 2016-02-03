%% produces supplementary figures for the journal version of our 
%  'criticality and correlations' project

clear all

% each figure is its own cell and can be run independently of the others, 
% given that the figure formatting variables defined in the following lines
% are in the workspace. 

axesThickness  = 1.0; % overall thickness of axes. Doesn't seem to do much.

fontName = 'Arial';    fontWeight     = 'normal';
fontSize       = 1 * 10;   fontSizeTitle  = 1 * 16;   
fontSizeXlabel = 1 * 10;   fontSizeYlabel = 1 * 11;
fontSizeText   = 1 * 10;   fontSizeLegend = 1 * 11;

%% supplementary figure 1: Advantage of Rao-Blackwellising 
figureS1 = figure('Tag', 'figS5', 'units','centimeters','position',...
                  [0,0,19,11]);

clrs = copper(20);
clrs = clrs(end:-1:1,:);

n = 100; % example population size chosen for the paper. 


for idxGroup = 1:3
  switch idxGroup
    case 1
      idxRange = 1:n;  
    case 2
      idxRange = (n+1) : ( n*(n+1)/2 ); % second moments / covariances        
    case 3
      idxRange = n*(n+3)/2+1 + (-n:0);  
  end
  timeRange = 1:65;     % select initial data segment for plotting, as the
                        % MSEs of the MCMC estimates necessarily become
                        % the same towards the end (we're comparing against
                        % the average of the final estimates from boths 
                        % MCMC chains)
  traces = cell(10, 2);                         
  finals = zeros( n * (n+3) / 2 +1, 10, 2);
  MSEperc = zeros(10,3);

 for idxRepet = 1:10
  tmp = cell(2,1);
  for r = 1:2
   switch r
    case 1
      load(['../results/method_validation/feature_moments_test_rao_blackwell/checkRB_RB100idxRepet', ...
             num2str(idxRepet), 'T1.mat'])
    case 2
      load(['../results/method_validation/feature_moments_test_rao_blackwell/checkRB_noRB100idxRepet', ...
             num2str(idxRepet), 'T1.mat'])
   end
   tmp{r} = Efy(:,2:end-1);
  
   fDescrJ = nchoosek(1:n,2)'; 
   covsy = tmp{r}((n+1):(n*(n+1)/2),:) - ...    % covariance as computed
          (tmp{r}(fDescrJ(1, :),:).* tmp{r}(fDescrJ(2, :),:)); % from Efx
   % now gives [firing rates; covariances; P(K)] over time 
   tmp{r}(n+(1:n*(n-1)/2),:) = covsy;  
   % moving average over time (i.e. tmp(:,t) holds MCMC estimate at time t)
   tmp{r} = bsxfun(@rdivide, cumsum(tmp{r},2), (1:size(tmp{r},2))); 
   % chain ends (for comparison with 'truth')
   finals(:, idxRepet, r) = tmp{r}(:,end);  
  end
  MSEperc(idxRepet,1) = mean( (finals(1:n, idxRepet, 1) ...
                             - finals(1:n, idxRepet, 2)).^2  ) / ...
                 mean( mean(finals(1:n, idxRepet, 1).^2) ...
                     + mean(finals(1:n, idxRepet, 2).^2));
  MSEperc(idxRepet,2) = mean( (finals(n+(1:n*(n-1)/2),idxRepet,1) ...
                             - finals(n+(1:n*(n-1)/2),idxRepet,2)).^2)/ ...
                 mean( mean(finals(n+(1:n*(n-1)/2), idxRepet, 1).^2) ...
                     + mean(finals(n+(1:n*(n-1)/2), idxRepet, 2).^2));
  MSEperc(idxRepet,3) = mean( (finals(end-n:end, idxRepet, 1) ...
                            - finals(end-n:end, idxRepet, 2)).^2  ) / ...
                 mean( mean(finals(end-n:end, idxRepet, 1).^2) ...
                     + mean(finals(end-n:end, idxRepet, 2).^2));
  for r = 1:2
   % compare to average of both estimates   
   tmp{r} = bsxfun(@minus, tmp{r}(:,:), mean(finals(:, idxRepet, :),3)).^2; 
   % compute average squared distance over all variables of interest
   traces{idxRepet, r} = (mean(tmp{r}(idxRange,:),1)' ...
                       ./ mean(mean(finals(idxRange,idxRepet,:),3).^2)); 
  end 
  clear tmp
 end

 avgtraces = cell(2,1);
 subplot(3,2,2*(idxGroup-1)+1)
 for r = 1:2
   for idxRepet = 1:10
    plot(100 * traces{idxRepet,r}(timeRange), ...
         'color', clrs(5*r,:), 'linewidth', 1.5);    
    hold on
   end
   
   m = Inf;
   for idxRepet = 1:10
       m = min([m, length(traces{idxRepet,r})]);
   end
   avgtraces{r} = zeros(m,1);
   for idxRepet = 1:10
     avgtraces{r} = avgtraces{r} + traces{idxRepet,r}(1:m);
   end 
   avgtraces{r} = avgtraces{r} / idxRepet;
 end
 for i = 1:5
  line([8, 8]*2^(i-1),[0, 1000],'linestyle','--','color',(i-1)/10*[1,1,1])
  hold on
 end
 hold off
 switch idxGroup
     case 1
       axis([0, 1.0*length(timeRange), 0, 12]) 
     case 2
       axis([0, 1.0*length(timeRange), 0, 200]) 
       ylabel('% normalised MSE in covariance from final estimate', ...
              'FontName', fontName, 'FontSize', fontSizeXlabel, ...
              'FontWeight', fontWeight ) 
     case 3
       axis([0, 1.0*length(timeRange), 0, 1]) 
       xlabel('# of sweeps', 'FontName', fontName, ...
             'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
 end
 set(gca, 'Linewidth', axesThickness)
 set(gca, 'TickDir', 'out'), box off
 set(gca, 'XTick', [8,16,32,64])
 set(gca, 'XTickLabel', 100*[8,16,32,64])
 set(gca, 'FontSize', fontSize)

 idxShow = unique(round(exp(-2:0.1:10)))+1; 
 idxShow = idxShow(idxShow<=length(avgtraces{1}));

 subplot(3,2,2*(idxGroup-1)+2)
 loglog(idxShow, 100 * avgtraces{1}(idxShow), 'color', clrs(5,:), 'linewidth',  2.5)
 hold on
 loglog(idxShow, 100 * avgtraces{2}(idxShow), 'color', clrs(10,:), 'linewidth', 2.5)
 loglog(idxShow, 100 * avgtraces{1}(idxShow), 'color', clrs(5,:), 'linewidth',  2.5)
 hold on
 loglog(idxShow, 100 * avgtraces{2}(idxShow), 'color', clrs(10,:), 'linewidth', 2.5)
 
 hold off
 switch idxGroup
     case 1
       axis([0.9, 0.525*m, 2*10^(-4), 5]) 
       set(gca, 'YTick', 10.^[-3, -2, -1, 0])
     case 2
       axis([0.9, 0.525*m, 9*10^(-4), 2.6*10^2]) 
       set(gca, 'YTick', 10.^[-2, 0, 2])
       ylabel('% normalised MSE rel. to final estimate', ...
              'FontName', fontName, 'FontSize', fontSizeXlabel, ...
              'FontWeight', fontWeight ) 
       legend('Rao-Blackwell', 'no Rao-Blackwell') 
     case 3
       axis([0.9, 0.525*m, 3*10^(-6), 0.2]) 
       set(gca, 'YTick', 10.^[-5, -3, -1])
       xlabel('# of sweeps', 'FontName', fontName, ...
              'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
 end
 set(gca, 'Linewidth', axesThickness)
 set(gca, 'TickDir', 'out'), box off
 set(gca, 'XTick', [2, 10, 50, 250, 1250, 6250])
 set(gca, 'XTickLabel', 100*[2, 10, 50, 250, 1250, 6250])
 set(gca, 'FontSize', fontSize)

end

%clearvars -except figureS1 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight

%% Compute gain of Rao-Blackwell over basic pair-wise Gibbs sampler
% parts of this subplot actually moved into the submission for the paper! 

figure42 = figure('Tag', 'fig42', 'units','centimeters','position',...
                    [0,0,19,11]);

clrs = copper(20);
clrs = clrs(end:-1:1,:);

n = 100; % example population size chosen for the paper. 

idxRepets = 1:10;

for idxGroup = 2
%  subplot(1,3,idxGroup)
  switch idxGroup
    case 1
      idxRange = 1:n;  
    case 2
      idxRange = (n+1) : ( n*(n+1)/2 ); % second moments / covariances        
    case 3
      idxRange = n*(n+3)/2+1 + (-n:0);  
  end
  timeRange = 1:65;                % first 15.000 sweeps
  traces = cell(10, 2);
  finals = zeros( n * (n+3) / 2 +1, 10, 2);
  MSEperc = zeros(10,3);

 for idxRepet = idxRepets
  tmp = cell(2,1);
  for r = 1:2
   switch r
    case 1
     load(['../results/method_validation/feature_moments_test_rao_blackwell/', ...
           'checkRB_RB100idxRepet', num2str(idxRepet), 'T1.mat'])
    case 2
     load(['../results/method_validation/feature_moments_test_rao_blackwell/', ...
           'checkRB_noRB100idxRepet', num2str(idxRepet), 'T1.mat'])
   end
   tmp{r} = Efy(:,2:end-1);
  
   fDescrJ = nchoosek(1:n,2)'; 
   covsy = tmp{r}((n+1):(n*(n+1)/2),:) - ...  % covariance as computed
          (tmp{r}(fDescrJ(1, :),:) ...        % from Efx    
        .* tmp{r}(fDescrJ(2, :),:)); 
   % now gives [firing rates; covariances; P(K)] over time 
   tmp{r}(n+(1:n*(n-1)/2),:) = covsy;  
   % moving average over time (i.e. tmp(:,t) holds MCMC estimate at time t)
   tmp{r} = bsxfun(@rdivide, cumsum(tmp{r},2), (1:size(tmp{r},2))); 
   % chain ends (for comparison with 'truth')
   finals(:, idxRepet, r) = tmp{r}(:,end);  
  end
  MSEperc(idxRepet,1) = mean( (finals(1:n, idxRepet, 1) ...
                             - finals(1:n, idxRepet, 2)).^2  ) / ...
                 mean( mean(finals(1:n, idxRepet, 1).^2) ...
                     + mean(finals(1:n, idxRepet, 2).^2));
  MSEperc(idxRepet,2) = mean( (finals(n+(1:n*(n-1)/2),idxRepet,1) ...
                             - finals(n+(1:n*(n-1)/2),idxRepet,2)).^2)/ ...
                 mean( mean(finals(n+(1:n*(n-1)/2), idxRepet, 1).^2) ...
                     + mean(finals(n+(1:n*(n-1)/2), idxRepet, 2).^2));
  MSEperc(idxRepet,3) = mean( (finals(end-n:end, idxRepet, 1) ...
                             - finals(end-n:end, idxRepet, 2)).^2  ) / ...
                 mean( mean(finals(end-n:end, idxRepet, 1).^2) ...
                     + mean(finals(end-n:end, idxRepet, 2).^2));
  for r = 1:2
   % compare to average of both estimates   
   tmp{r} = bsxfun(@minus, tmp{r}(:,:), mean(finals(:, idxRepet, :),3)).^2; 
   % compute average squared distance over all variables of interest
   traces{idxRepet, r} = (mean(tmp{r}(idxRange,:),1)' ...
                         ./ mean(mean(finals(idxRange,idxRepet,:),3).^2)); 
  end 
 end

 % compute averages over all 10 runs each for RB and no RB condition
 avgtraces = cell(2,1);
 stdtraces = cell(2,1);
 for r = 1:2
   m = Inf;
   for idxRepet = 1
       m = min([m, length(traces{idxRepet,r})]);
   end
   avgtraces{r} = zeros(m,1);
   stdtraces{r} = zeros(m,length(idxRepets));
   for idxRepet = idxRepets
     avgtraces{r} = avgtraces{r} + traces{idxRepet,r}(1:m);        % MSE
     stdtraces{r}(:,idxRepet) = traces{idxRepet,r}(1:m);           % MSE
%     avgtraces{r} = avgtraces{r} + sqrt(traces{idxRepet,r}(1:m)); % RMSE
%     stdtraces{r}(:,idxRepet) = sqrt(traces{idxRepet,r}(1:m));    % RMSE
   end 
   avgtraces{r} = avgtraces{r} / idxRepet;
   stdtraces{r} = std(stdtraces{r},0,2);
 end
 
idxShow = unique(round(exp(-2:0.1:10)))+1; 
idxShow = idxShow(idxShow<=length(avgtraces{1}));

loglog([0.00001,0.00001], '-', 'color', 'k', 'linewidth',  2.5)
hold on
loglog([0.00001,0.00001], ':', 'color', 'k', 'linewidth',  2.5)

h = area(idxShow, ...
     [ 100 * avgtraces{2}(idxShow) - 100 * stdtraces{2}(idxShow), ...
       200 * stdtraces{2}(idxShow)]);
h(1).FaceColor = [1,1,1];
h(1).EdgeColor = 'none';
h(2).FaceColor = [clrs(7-2*idxGroup,:)];
h(2).EdgeColor = 'none';
loglog(idxShow, 100 * avgtraces{2}(idxShow), ...
       ':', 'color', 'k', 'linewidth',  2.5)

h = area(idxShow, ...
     [ 100 * avgtraces{1}(idxShow) - 100 * stdtraces{1}(idxShow), ...
       200 * stdtraces{1}(idxShow)]);
h(1).FaceColor = [1,1,1];
h(1).EdgeColor = 'none';
h(2).FaceColor = [clrs(7-2*idxGroup,:)];
h(2).EdgeColor = 'none';
loglog(idxShow, 100 * avgtraces{1}(idxShow), ...
       '-', 'color', 'k', 'linewidth',  2.5)

switch idxGroup
    case 1
        axis([0, m/2, 0.08, 20])
    case 2
        axis([0, m/2, 0.01, 320])
    case 3
        axis([0, m/2, 0.009, 3])
end
set(gca, 'Linewidth', axesThickness)
set(gca, 'TickDir', 'out'), box off
set(gca, 'XTick', [2, 10, 50, 250, 1250, 6250])
set(gca, 'XTickLabel', 100*[2, 10, 50, 250, 1250, 6250])
set(gca, 'FontSize', fontSize)
if idxGroup == 1
 ylabel('% normalised RMSE rel. to final estimate', ...
       'FontName', fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight ) 
end
xlabel('# of sweeps', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight ) 
       legend('Rao-Blackwell', 'no Rao-Blackwell') 

end
%% Quantify gain of Rao-Blackwellising
%  (figure not to appear in the paper)

% we do two pairs of linear approximations in the log-log space: 
% 1) two lines with invidivual slopes (log(y) = a1/2 * log(x) + b1/2
% 2) two lines with same slopes (log(y) = a3 * log(x) + b3_1/2

x = 100*(8:300); y1 = avgtraces{1}(8:300)'; y2 = avgtraces{2}(8:300)';
logy1 = log(y1);  logy2 = log(y2); logy3 = [logy2(:);logy1(:)];              
logx = log(x); logx3 = [logx(:); logx(:)];

% allow two different slopes a1, a2 for both lines
a1 = (sum(logy1)*sum(logx) - length(logy1) * sum(logx.*logy1)) ...
      / (sum(logx)^2 - length(logy1) *sum(logx.^2));
a2 = (sum(logy2)*sum(logx) - length(logy2) * sum(logx.*logy2)) ...
      / (sum(logx)^2 - length(logy2) *sum(logx.^2));
b1 = (sum(logy1) - a1 * sum(logx))/length(logy1);
b2 = (sum(logy2) - a2 * sum(logx))/length(logy2);
% enforce identical slope a3 on both lines
a3 = (sum(logy3)*sum(logx3) - length(logy3) * sum(logx3.*logy3)) ...
      / (sum(logx3)^2 - length(logy3) *sum(logx3.^2));
b31 = (sum(logy1) - a3 * sum(logx))/length(logy1);
b32 = (sum(logy2) - a3 * sum(logx))/length(logy2);

figure(3919); 
% plot the averaged MSE for both conditions, and the linear approximations
% on top of it (ugly colors)
subplot(1,3,1)
        plot(log(100*(1:length(avgtraces{1}))),log(avgtraces{1}), ...
              'k','linewidth',2); 
        hold on; plot(logx, a1*logx+b1, 'go-','linewidth',2); hold on
        plot(log(100*(1:length(avgtraces{2}))),log(avgtraces{2}), ...
            'k','linewidth',2); 
        hold on; plot(logx, a2*logx+b2, 'ro-','linewidth',2);
        hold on; plot(logx, a3*logx+b31, 'mx-','linewidth',2);
        hold on; plot(logx, a3*logx+b32, 'cx-','linewidth',2);
        xlabel('log(#sweeps)')
        ylabel('log(MSE)')
        
subplot(1,3,2); 
% plot the ratio of normalised MSE as a function of #sweeps, and the
% results from the linear approximations on
plot(x,avgtraces{2}(8:300)'./avgtraces{1}(8:300)'); 
hold on; 
plot(x, ones(size(logx))*exp(b32-b31), 'c'); 
plot(x, exp(b2-b1) * (x).^(a2-a1),'r')     
xlabel('#sweeps')
ylabel('MSE_{no RB}/MSE_{RB}')
legend('data', 'linear fits, a_1=a_2', 'linear fits with a_1, a_2')

subplot(1,3,3); 
% plot the ratio of needes samples for the same MSE as a function of 
% #sweeps both for RB and no RB condition. Also overlay the results from
% the line fits. 
closest1 = zeros(size(y1));
check = zeros(size(y1));
for i = 1:length(y1)
    [~,idx] = min(abs(y1- y2(i))); % for each RB ...
    check(i) = y1(idx);
    closest1(i) = x(idx);        % find closest with no RB
end
plot(x, exp(logx)./closest1, 'b') % real data
hold on

t1 = a3 * logx + b31; % RB
t2 = a3 * logx + b32; % no RB
closest1 = zeros(size(t2));
for i = 1:length(t1)
    [~,idx] = min(abs(t1-t2(i))); % for each RB ...
    closest1(i) = logx(idx);        % find closest with no RB
end
plot(x, exp(logx)./exp(closest1), 'r')
plot(x, ones(size(logx))*exp((b31-b32)/a3), 'm');   

% that one is a bit of a mess: 
% I found it way easer to flip the ratios for this case, i.e. this one
% plots the number of RB sweeps needed to no-RB-performance. The important
% value is the ratio of around (roughly) 30% of samples needed. 
t1 = a1 * logx + b1; % RB
t2 = a2 * logx + b2; % no RB
closest1 = zeros(size(t2));
for i = 1:length(t1)
    [~,idx] = min(abs(t1- t2(i))); % for each RB ...
    closest1(i) = logx(idx);        % find closest with no RB
end
plot(x, exp(closest1)./x, 'k')
hold on; 
plot(x, exp((b2-b1)/a1) * x.^(a2/a1-1),'c')   
set(gca, 'XTick', 800*2.^[0:2:8])
xlabel('#sweeps')
ylabel('#samples_{no RB}/#samples_{RB}')
legend('data', 'closest data', 'linear fits, a_1=a_2', ...
               'closest data', 'linear fits with a_1, a_2')
           
clearvars -except figure42 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight           
           
%% supplementary figure 2: convergence and quality of K-pairwise model
%                          for 10 fits at population size n = 100

figureS1 = figure('Tag', 'figS1', 'units','centimeters','position',...
                  [0,0,19,11]);
n = 100;
load(['../results/method_validation/K_pairwise_fit_n100/nat_n', num2str(n), 'asfinal'])  
load(['../results/K_pairwise_final/Efx_nat.mat'])
mx = 0;

clrs = [254,153,41;236,112,20;204,76,2;153,52,4;102,37,6;0,0,0]/255;

subplot(21,34,vec(bsxfun(@plus,  (1:9)', (1:20)*34)))

for i = 1:length(fD), 
    mx = max([mx, 10*(length(fD{i}.MSEperc(1,2:end-1)))]);
end
tmp = zeros(4, mx/10);

% switch between plotting RMSEs and MSEs: 
% The 'MSEs' stored in fD actually are RMSEs by now (need to change naming
% at some point...). Switch back to plotting MSEs simply by taking the 
% square root (^k, k=2), or keep RMSEs by applying identity (^k, k =1)
k = 2;

for i = 1:length(fD), 
    tmp(1,1:size(fD{i}.MSEperc,2)-2) = tmp(1,1:size(fD{i}.MSEperc,2)-2) ...
               + fD{i}.MSEperc(1,2:end-1).^k;
    tmp(2,1:size(fD{i}.MSEperc,2)-2) = tmp(2,1:size(fD{i}.MSEperc,2)-2) ...
               + fD{i}.MSEperc(4,2:end-1).^k;
    tmp(3,1:size(fD{i}.MSEperc,2)-2) = tmp(3,1:size(fD{i}.MSEperc,2)-2) ...
               + fD{i}.MSEperc(3,2:end-1).^k;
    tmp(4,1:size(fD{i}.MSEperc,2)-2) = tmp(4,1:size(fD{i}.MSEperc,2)-2) +1;
    plot(10*(1:length(fD{i}.MSEperc(1,2:end-1))),...
         100*fD{i}.MSEperc(1,2:end-1).^k, '.', 'color', clrs(5,:)); 
    hold on     
    plot(10*(1:length(fD{i}.MSEperc(1,2:end-1))),...
         100*fD{i}.MSEperc(3,2:end-1).^k, '.', 'color', clrs(1,:)); 
    plot(10*(1:length(fD{i}.MSEperc(1,2:end-1))),...
         100*fD{i}.MSEperc(4,2:end-1).^k, '.', 'color', clrs(3,:)); 

end
tmp = bsxfun(@rdivide, tmp(1:3,:), tmp(4,:));
mvgavg = zeros(3, mx/10);
kernel_length = 15; % correction mask below only accurate for uneven ...
                    % kernel lengths...
mask = [ kernel_length./(ceil(kernel_length/2):kernel_length-1), ...
         ones(1, mx/10 - (kernel_length-1)), ...
         kernel_length./(kernel_length-1:-1:ceil(kernel_length/2)) ];        
     
% conv(u,v) uses zero-padding at the borders of u, which pushes down the
% estimates. We correct for that by a multiplicative mask. The result
% behaves as if the padding was done with the average of those values of u
% still within the valid range covered by v.

mvgavg(1,:) = conv(tmp(1,:), 1/kernel_length*ones(1,kernel_length),'same');
mvgavg(2,:) = conv(tmp(2,:), 1/kernel_length*ones(1,kernel_length),'same');
mvgavg(3,:) = conv(tmp(3,:), 1/kernel_length*ones(1,kernel_length),'same');
plot(10*(1:mx/10),100*mvgavg(1,:).*mask, '-', 'color', clrs(2,:), ...
    'linewidth', 4); hold on, 
plot(10*(1:mx/10),100*mvgavg(2,:).*mask, '-', 'color', clrs(6,:), ...
    'linewidth', 4); 
plot(10*(1:mx/10),100*mvgavg(3,:).*mask, '-', 'color', clrs(6,:), ...
    'linewidth', 4); 
plot(10*(1:mx/10),100*mvgavg(1,:).*mask, '-', 'color', clrs(5,:), ...
    'linewidth', 2);
plot(10*(1:mx/10),100*mvgavg(2,:).*mask, '-', 'color', clrs(3,:), ...
    'linewidth', 2); 
plot(10*(1:mx/10),100*mvgavg(3,:).*mask, '-', 'color', clrs(1,:), ...
    'linewidth', 2);

set(gca, 'Linewidth', axesThickness)
box off, set(gca, 'TickDir', 'out')
set(gca, 'FontSize', fontSize)
set(gca, 'XTick', 2000:2000:6000)
xlabel('# sweeps', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight )
ylabel('% normalised MSE', 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
axis([5, mx+5, 0, 120])

subplot(21,34,vec(bsxfun(@plus,  (13:18)', (2:8)*34)))
for i = 1:length(fD), 
  plot(50*Efx{n/10}(1:n, i), 50*fD{i}.Efy(1:n), '.', 'color', clrs(5,:)) 
  hold on;
end
box off, set(gca, 'TickDir', 'out')
set(gca, 'FontSize', fontSize)
xlabel('true FR [Hz]',  'FontName',fontName, 'FontSize', fontSizeXlabel,... 
       'FontWeight', fontWeight )
ylabel('model estimate', 'FontName',fontName,'FontSize',fontSizeXlabel, ...
       'FontWeight', fontWeight )
set(gca, 'XTick', [0, 5])
set(gca, 'YTick', [0, 5])
axis([0, 6, 0, 6])

subplot(21,34,vec(bsxfun(@plus,  (21:26)', (2:8)*34)))
fDescrJ = nchoosek(1:n,2)'; 
for i = 1:length(fD), 
    covsx = Efx{n/10}((n+1):(n*(n+1)/2),i) - ... % covariance as computed
           (Efx{n/10}(fDescrJ(1, :),i).* ...
            Efx{n/10}(fDescrJ(2, :),i)); % from Efx    for j = js, 
    varsx = Efx{n/10}(1:n,i) .* (1 - Efx{n/10}(1:n,i));
    corsx = covsx ./(sqrt(varsx(fDescrJ(1,:))).*sqrt(varsx(fDescrJ(2,:))));
    Efy = fD{i}.Efy;
    covsy = Efy((n+1):(n*(n+1)/2)) - ... % covariance as computed
           (Efy(fDescrJ(1, :)).* Efy(fDescrJ(2, :))); % from Efx    
    varsy = Efy(1:n) .* (1 - Efy(1:n));
    corsy = covsy ./(sqrt(varsy(fDescrJ(1,:))).*sqrt(varsy(fDescrJ(2,:))));
  plot(covsx, covsy, '.', 'color', clrs(3,:)), hold on;
end
box off, set(gca, 'TickDir', 'out')
set(gca, 'FontSize', fontSize)
xlabel('true covariances',  'FontName',fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )
set(gca, 'XTick', [0, 0.05])
set(gca, 'YTick', [0, 0.05])
axis([-0.01, 0.06, -0.01, 0.06])


subplot(21,34,vec(bsxfun(@plus,  (29:34)', (2:8)*34)))
for i = 1:length(fD), 
  plot(Efx{n/10}(end-n:end, i), fD{i}.Efy(end-n:end), '.', ...
       'color', clrs(1,:)), hold on;
end
box off, set(gca, 'TickDir', 'out')
set(gca, 'FontSize', fontSize)
xlabel('true P(K=k)',  'FontName',fontName, 'FontSize', fontSizeXlabel, ...
       'FontWeight', fontWeight )
set(gca, 'XTick', [0, 0.4])
set(gca, 'YTick', [0, 0.4])
axis([0, 0.48, 0, 0.48])

trues = {'cb', 'nat', 'fff'}; 
MSEs = cell(size(trues));
for j = 1:3
    
  load(['../results/method_validation/K_pairwise_fit_n100/', trues{j},'_n', num2str(n), 'asfinal'])  
  load(['../results/K_pairwise_final/Efx_', trues{j}, '.mat'])

  MSEs{j} = zeros(length(fD),3);
  for i = 1:length(fD)
   Efy = fD{i}.Efy;
   covsy = Efy((n+1):(n*(n+1)/2)) - ...            % covariance as computed
          (Efy(fDescrJ(1,:)).* Efy(fDescrJ(2,:))); % from Efx    
   Efy(n+(1:n*(n-1)/2)) = covsy;   
   tmp = Efx{10}(:,i);
   covsx = tmp((n+1):(n*(n+1)/2)) - ...            % covariance as computed
          (tmp(fDescrJ(1,:)).* tmp(fDescrJ(2,:))); % from Efx    
   tmp(n+(1:n*(n-1)/2)) = covsx;         
   for idxGroup = 1:3
    switch idxGroup
     case 1
       idxRange = 1:n;  
     case 2
       idxRange = (n+1) : ( n*(n+1)/2 ); % second moments / covariances        
     case 3
       idxRange = n*(n+3)/2+1 + (-n:0);  
     end
    %MSEs{j}(i,idxGroup)=100*sqrt(mean((Efy(idxRange)-tmp(idxRange)).^2)...
    %                    /mean(tmp(idxRange).^2)); % RMSE
    MSEs{j}(i,idxGroup)=100*(mean((Efy(idxRange)-tmp(idxRange)).^2) ...
                        /mean(tmp(idxRange).^2)); % MSE
   end
  end
end

subplot(21,34,vec(bsxfun(@plus,  (13:34)', (12:20)*34)))

for j = 1:3
plot(j-0.15, MSEs{j}(:,1), '.', 'color', clrs(5,:)), hold on;    
plot(j     , MSEs{j}(:,2), '.', 'color', clrs(3,:));    
plot(j+0.15, MSEs{j}(:,3), '.', 'color', clrs(1,:));    
line([j-0.05, j+0.05]-0.15, mean(MSEs{j}(:,1)) * [1,1], ...
     'color', clrs(5,:), 'linewidth', 1)    
line([j-0.05, j+0.05]     , mean(MSEs{j}(:,2)) * [1,1], ...
     'color', clrs(3,:), 'linewidth', 1)    
line([j-0.05, j+0.05]+0.15, mean(MSEs{j}(:,3)) * [1,1], ...
     'color', clrs(1,:), 'linewidth', 1)    
end
box off, set(gca, 'TickDir', 'out')
axis([0.5, 3.5, 0, 6])
xlabel('stimulus condition')
ylabel('% normalised MSE')
set(gca, 'XTick', 1:3)
set(gca, 'XTickLabel', {'cb', 'nat', 'fff'})

clearvars -except figureS1 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight

%% supplementary figure 3: randomized P(K) and retained criticality
figureS3 = figure('Tag', 'figS3', 'units','centimeters','position', ...
                  [0,0,19,11]);
load('fig_data/figS3_data_shuffle.mat')

clrs = copper(18);

% plot P(K) for full N = 300 population
subplot(2,3,1:2),  
idxUnshuffle = zeros(size(idxShuffle));
for i = 1:length(pcount)
  idxUnshuffle(i) = find(idxShuffle==i);
end
plot(0:length(pcount)-1, pcount(idxUnshuffle), 'color', clrs(10,:), ...
     'linewidth', 2.5)
hold on
plot(0:length(pcount)-1, pcount, 'k', 'linewidth', 2.5)
axis([0,length(pcount)-1, 0, 1.05*max(pcount)]);
set(gca, 'TickDir', 'out'), box off
set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick', [0, 100,200,300]);
set(gca, 'YTick', [0, 0.1, 0.2])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('population spike count'   , 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('probability', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
       'FontWeight', fontWeight )

% add what the P(K) actually look like for randomly subsampled n < 300
subplot(2,3,4:5),
for i = 1:15, 
    plot(0:size(pcounts,1)-1, squeeze(pcounts(:,i,:)), ...
         'color', clrs(19-i,:)); 
    hold on, 
end
axis([0,length(pcount)-1, 0, 1.05*max(pcounts(:))]);
set(gca, 'TickDir', 'out'), box off
set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick', [0, 100,200,300]);
set(gca, 'YTick', [0, 0.04, 0.08])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('population spike count'   , 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('probability', 'FontName', fontName, 'FontSize', fontSizeYlabel, ...
       'FontWeight', fontWeight )

% add heat traces computed from above subsampled populations
subplot(2,3,[3,6]), 
lgnd = cell(15,1);
for i = 15:-1:1
    plot(-1, 0, 'color', clrs(19-i, :), 'linewidth', 1.5); 
    hold on
    lgnd{i} = ['n = ', num2str(Ns(16-i))];
end
for i = 1:15, 
    plot(Ts, squeeze(cN(i,:,:)), 'color', clrs(19-i,:)); 
    hold on, 
end
line([1,1], [0, 1.05*max(cN(:))], 'linestyle', '--', 'color', 'k', ...
     'linewidth', axesThickness)
axis([min(Ts), max(Ts), 0.95*min(cN(:)), 1.05*max(cN(:))]);
set(gca, 'TickDir', 'out'), box off
set(gca, 'Linewidth', axesThickness)
set(gca, 'XTick', [0.7, 1, 1.2]);
set(gca, 'YTick', [0, 5, 10])
box off, set(gca, 'TickDir' ,'out')
legend(lgnd, 'Location', 'Northwest'), legend boxoff
set(gca, 'FontSize', fontSize)
xlabel('temperature'   , 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('specific heat', 'FontName', fontName, ...
       'FontSize', fontSizeYlabel, 'FontWeight', fontWeight )

clearvars -except figureS3 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight

%% supplementary figure 4: c(T) traces for small simulation up to n = 300, with and without beta-bin. approx.

figureS4 = figure('Tag', 'figS4', 'units','centimeters','position',[0,0,19,11]);

lgnd = cell(15,1);
clrs = copper(15);
clrs = clrs(end:-1:1,:);

subplot(1,2,1)
load('fig_data/figS4_data.mat')
for i = 15:-1:1, 
    plot(Ts, mean(squeeze(cN(i,:,:)),1), 'color', clrs(i,:), ...
         'linewidth', 2.5), hold on; 
    lgnd{16-i} = ['n = ', num2str(Ns(i))];
end
for i = 15:-1:1, 
    plot(Ts, squeeze(cN(i,:,:))', 'linestyle', '-', 'linewidth', 1, ...
         'color', clrs(i,:)), 
end

line([1,1], [0, 1.05*max(cN(:))], 'linestyle', '--', 'color', 'k', ...
     'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([0.95,1.2,0.95*min(cN(:)), 1.05*max(cN(:))]); %axis autoy
set(gca, 'XTick', [1, 1.1, 1.2, 2]);
set(gca, 'YTick', [5,10])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
xlabel('temperature'   , 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('specific heat', 'FontName', fontName, ...
       'FontSize', fontSizeYlabel, 'FontWeight', fontWeight )

subplot(1,2,2)
clear cN Ns Ts
load('fig_data/figS4_data_raw.mat')
for i = 15:-1:1, 
    plot(Ts, mean(squeeze(cN(i,:,:)),1), 'color', clrs(i,:), ...
         'linewidth', 2.5), hold on; 
    lgnd{16-i} = ['n = ', num2str(Ns(i))];
end
for i = 15:-1:1, 
    plot(Ts, squeeze(cN(i,:,:))', 'linestyle', '-', 'linewidth', 1, ...
         'color', clrs(i,:)), 
end

line([1,1], [0, 1.05*max(cN(:))], 'linestyle', '--', 'color', 'k', ...
     'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([0.95,1.2,0.95*min(cN(:)), 1.05*max(cN(:))]); %axis autoy
set(gca, 'XTick', [1, 1.1, 1.2, 2]);
set(gca, 'YTick', [5,10])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
legend(lgnd); legend boxoff
xlabel('temperature'   , 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('specific heat', 'FontName', fontName, ...
       'FontSize', fontSizeYlabel, 'FontWeight', fontWeight )

clearvars -except figureS4 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight

%% Key population statistics (FR, covs, P(K)) as function of temperature
figureS9 = figure('Tag', 'figS8', 'units','centimeters','position',[0,0,19,11]);
Ts =   [0.8000,    0.8563,    0.9000,    0.9363,    0.9675,    0.9938, ...   
        1.0175,    1.0388,    1.0575,    1.0775,    1.0988,    1.1200, ...
        1.1413,    1.1650,    1.1900,    1.2150,    1.2425,    1.2713, ...   
        1.3013,    1.3338,    1.3687,    1.4063,    1.4463,    1.4900, ...
        1.5375,    1.5900,    1.6500,    1.7175,    1.7950,    1.8875, ...   
        2.0000,    1.0000];
[~, idx] = sort(Ts);

Tis = [1:5, 7:32];
clrs = jet(31);
n = 100;
fDescrJ = nchoosek(1:n,2)';  

Temps = zeros(31,1);
for i = 1:31, 
    Temps(i) = Ts(idx(Tis(i)));
end
for i = 1:31
    str = num2str(Temps(i));
    disp(str)
    str(str=='.') = '_';
    
%    load(['../results/method_validation/specific_heat_samples_long_runs/','longRun100idxRepet1T',str,'.mat']) 
    load(['/home/mackelab/Desktop/Projects/Criticality/results/VarE_long_runs/longRun100idxRepet1T',str,'.mat'])
    
    subplot(1,3,1)
    plot(Temps(i),mean(50*Efy(1:n)),'.','color',clrs(i,:),'markerSize',10)
    hold on    
    line(Temps(i)*[1,1],mean(50*Efy(1:n))+std(50*Efy(1:n))/ ...
                                                       sqrt(n)*[-1,1], ...
         'color', clrs(i,:))
    set(gca, 'TickDir', 'out'), box off
    axis([0.75, 2.05, 0, 1]), 
    axis autoy
    set(gca, 'XTick', [1, 1.5, 2.0])
    set(gca, 'YTick', [0, 5, 10, 15, 20])
    set(gca, 'Linewidth', axesThickness)
    set(gca, 'FontSize', fontSize)        
    xlabel('temperature', 'FontName', fontName, ...
           'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )  
    ylabel('FR [Hz]', 'FontName', fontName, 'FontSize', fontSizeXlabel, ...
           'FontWeight', fontWeight )  

    subplot(1,3,2)
    covsy = Efy((n+1):(n*(n+1)/2)) - ...       % covariance as computed
            (Efy(fDescrJ(1, :)).* Efy(fDescrJ(2, :))); % from Efx     
    vary  = Efy(1:n) .* (1 - Efy(1:n));

    plot(Temps(i), mean(covsy), '.', 'color', clrs(i,:), 'markerSize', 10)
    line(Temps(i)*[1,1], mean(covsy)+std(covsy)/ ...
                            sqrt(length(covsy))*[-1,1], 'color', clrs(i,:))
    hold on    
    set(gca, 'TickDir', 'out'), box off
    axis([0.75, 2.05, 0, 1]), 
    axis autoy
    set(gca, 'XTick', [1, 1.5, 2.0])
    set(gca, 'YTick', [0, 0.03, 0.06, 0.09])
    set(gca, 'YTick', [-0.01, 0, 0.01])
    set(gca, 'Linewidth', axesThickness)
    set(gca, 'FontSize', fontSize)        
    xlabel('temperature', 'FontName', fontName, ...
           'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )  
    ylabel('covariances', 'FontName', fontName, ...
           'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )  

    subplot(1,3,3)
    if Temps(i) == 1
        semilogy(0:n, PK(:, end-1), 'color', [254,153,41]/256, ...
                 'linewidth', 3.5)
        hold on
    end
    semilogy(0:n, PK(:, end-1), 'color', clrs(i,:))
    hold on    
    set(gca, 'TickDir', 'out'), box off
    axis([0, 47, 5*10^(-7), 1])    
    set(gca, 'XTick', [0, 20, 40])
    %set(gca, 'YTick', 10.^[-6, -4, -2, 0])
    set(gca, 'Linewidth', axesThickness)
    set(gca, 'FontSize', fontSize)        
    xlabel('population spike count', 'FontName', fontName, ...
           'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )  
    ylabel('probability', 'FontName', fontName, ...
           'FontSize', fontSizeXlabel, 'FontWeight', fontWeight )  

end

clearvars -except figureS9 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight

%% c(T) traces for large simulation, up to n = 6000, uniform vs. spatial subsampling

figureS3 = figure('Tag', 'figS3', 'units','centimeters','position',[0,0,19,11]);
lgnd = cell(30,1);
clrs = copper(30);
clrs = clrs(end:-1:1,:);

clear Ts Ns cN
load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Figures/Fig4/fig4rnd_data.mat')
cNrnd = cN; clear cN;
load('/home/mackelab/Desktop/Projects/Criticality/code_critical_retina/Figures/Fig4/fig4lin_data.mat')
cNlin = cN; clear cN;

offset = 0.2;
tempRange = 51:190;
for i = 30:-1:1, 
    plot(Ts(tempRange), mean(squeeze(cNrnd(10*i,:,tempRange)),1), ...
         'color', clrs(i,:), 'linewidth', 1), hold on; 
    lgnd{31-i} = ['n = ', num2str(Ns(10*i))];
end
for i = 30:-1:1, 
    plot(Ts(tempRange)+offset, squeeze(cNlin(10*i,:,tempRange)), ...
         'color', clrs(i,:), 'linewidth', 1), hold on; 
end
for i = 30:-1:1,
    [maxCrnd, idxTmaxCrnd] = max( mean(squeeze(cNrnd(10*i,:,:)),1) );
    [maxClin, idxTmaxClin] = max( squeeze(cNlin(10*i,:,:)) );
end
legend(lgnd); 
line([1,1], [0, 1.05*max(cNlin(:))], 'linestyle', '--', 'color', 'k', ...
     'linewidth', axesThickness)
line([1,1]+offset, [0, 1.05*max(cNlin(:))], 'linestyle', '--', ...
     'color', 'k', 'linewidth', axesThickness)
set(gca, 'Linewidth', axesThickness)
axis([min(Ts),max(Ts),0.95*min(cNlin(:)), 1.05*max(cNlin(:))]); %axis autoy
set(gca, 'XTick', [1, 1.1, 1.2, 1.3]);
set(gca, 'XTickLabel', [1, 1.1, 1, 1.1]);
set(gca, 'YTick', [10,20,30])
box off, set(gca, 'TickDir' ,'out')
set(gca, 'FontSize', fontSize)
legend(lgnd, 'orientation', 'vertical'); %legend boxoff
xlabel('temperature'   , 'FontName', fontName, ...
       'FontSize', fontSizeXlabel, 'FontWeight', fontWeight ) 
ylabel('specific heat', 'FontName', fontName, ...
       'FontSize', fontSizeYlabel, 'FontWeight', fontWeight )

clearvars -except figureS3 axesThickness clrs fontName fontSize fontSizeLegend fontSizeText fontSizeTitle fontSizeXlabel fontSizeYlabel fontWeight

