function [Ef,Es] = testFormulasBetaDistr(as, bs, ms, ns, ks, ls, nSamples, ifplot)

if nargin < 8
  ifplot = false;
end
if nargin < 7
  nSamples = 100000; 
end 

if nargin < 1
  as = [0.05:0.05:2];
end
if nargin < 2
  bs = [0.05:0.05:2];
end
if nargin < 3
  ms = [2]; % exponents m of log(x)^m
end
if nargin < 4
  ns = [0]; % exponents n of log(1-x)^n
end
if nargin < 5
  ks = [2]; % exponents k of x^k
end
if nargin < 6
  ls = [0]; % exponents l of x^l
end


Ef = zeros(length(ms),length(ns),length(ks),length(ls),length(as),length(bs));
Es = zeros(size(Ef));
for idxa = 1:length(as)
 a = as(idxa);    
 for idxb = 1:length(bs)
  b = bs(idxb);      
  x = betarnd(a,b,nSamples,1);
  for idxm = 1:length(ms)
   m = ms(idxm);       
   for idxn = 1:length(ns)
    n = ns(idxn);    
    for idxk = 1:length(ks)
     k = ks(idxk);    
     for idxl = 1:length(ls)
      l = ls(idxl);
      Ef(idxm,idxn,idxk,idxl,idxa,idxb) = prefactor(k,l,a,b) * g(m,n,a+k,b+l);    
      Es(idxm,idxn,idxk,idxl,idxa,idxb) = mean(log(x).^m .* x.^k .* log(1-x).^n .* (1-x).^l);
     end
    end
   end
  end
 end
end


function [out] = prefactor(k,l,a,b) 
  out = prod(a+(0:(k-1))) * prod(b+(0:(l-1))) / prod(a+b+(0:(k+l-1)));
end

 function [out] = g(m,n,a,b) 
  if m==1 && n==0
    out = psi(0,a) - psi(0,a+b);  
  elseif m==0 && n==1
    out = psi(0,b) - psi(0,a+b);  
  elseif m==1 && n==1
    out = psi(0,a)*psi(0,b) - psi(0,a+b)*(psi(0,a)+psi(0,b)) ...
        + psi(0,a+b)^2 - psi(1,a+b);
  elseif m==0 && n==2
    out = (psi(0,b) - psi(0,a+b))^2 + psi(1,b) - psi(1,a+b);
  elseif m==2 && n==0
    out = (psi(0,a) - psi(0,a+b))^2 + psi(1,a) - psi(1,a+b);
     
  else
    out = 1;
    warning('Warning: g_{(m,n)}(\alpha,\beta) not implemented for current (m,n)!')
  end

 end

 if ifplot
  figure; 
  subplot(1,2,1); 
  imagesc(squeeze(Ef)); 
  xlabel('\alpha');ylabel('\beta')
  title(['Analytic results for E[log(x)^', num2str(ms),' x^', num2str(ks), ...
         ' log(1-x)^', num2str(ns),'  (1-x)^', num2str(ls),']'])
  set(gca, 'XTick', round(linspace(1,length(as), 4)));
  set(gca, 'XTickLabel', as(round(linspace(1,length(as), 4))));
  set(gca, 'YTick', round(linspace(1,length(bs), 4)));
  set(gca, 'YTickLabel', bs(round(linspace(1,length(bs), 4))));
  set(gca, 'TickDir', 'out'); box off
  subplot(1,2,2); 
  imagesc(squeeze(Es));
  xlabel('\alpha');ylabel('\beta')
  title(['Monte Carlo results for E[log(x)^', num2str(ms),' x^', num2str(ks), ...
         ' log(1-x)^', num2str(ns),'  (1-x)^', num2str(ls),']'])
  set(gca, 'XTick', round(linspace(1,length(as), 4)));
  set(gca, 'XTickLabel', as(round(linspace(1,length(as), 4))));
  set(gca, 'YTick', round(linspace(1,length(bs), 4)));
  set(gca, 'YTickLabel', bs(round(linspace(1,length(bs), 4))));
  set(gca, 'TickDir', 'out'); box off

  maxRelDiff=max(vec(abs(squeeze(Ef)-squeeze(Es)))./abs(vec(squeeze(Ef))));
  annotation('textbox', [0 0.9 1 0.1], ...
    'String', ['Maximal relative error max\{ (E_{pred} - E_{MC}) / E_{pred} \} =  ', num2str(round(maxRelDiff*1000)/10), '%'], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')  
 end
end