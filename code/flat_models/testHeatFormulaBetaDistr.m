function [cf,cs,ci] = testHeatFormulaBetaDistr(as, bs, nSamples, ifplot)

if nargin < 4
  ifplot = false;
end
if nargin < 3
  nSamples = 100000; 
end 

if nargin < 1
  as = [0.05:0.05:2];
end
if nargin < 2
  bs = [0.05:0.05:2];
end

cf = zeros(length(as),length(bs));
cs = zeros(size(cf));
ci = zeros(size(cf));
for idxa = 1:length(as)
 a = as(idxa);    
 for idxb = 1:length(bs)
  b = bs(idxb);      
  x = betarnd(a,b,nSamples,1);
  %cf(idxa, idxb) = h(a,b);
  %cs(idxa, idxb) = mean(eta(x));
  %cf(idxa,idxb) =   prefactor(2,0,a,b)*g(2,0,a+2,b  ) ...
  %              + 2*prefactor(1,1,a,b)*g(1,1,a+1,b+1) ...
  %              +   prefactor(0,2,a,b)*g(0,2,a,  b+2) ...
  %              - h(a,b)^2;
  cf(idxa,idxb) = ( a*(a+1)*psi(1,a+1) + b*(b+1)*psi(1,b+1) )/((a+b)*(a+b+1)) ...
                - psi(1,a+b+1) ...
                + a*b*(psi(0,a+1)-psi(0,b+1))^2/((a+b)^2*(a+b+1));
  cs(idxa,idxb)  = mean( (eta(x) - h(a,b)).^2 ); 

  fh = @(x) betapdf(x,a,b) .* (eta(x) - h(a,b)).^2; % function handle
  dx = 0.000001;
  ip = 0.5*dx:dx:1;            % integration points
  ci(idxa,idxb) = sum(fh(ip))*dx;
%  ci(idxa,idxb) = integral(fh, 0, 1);
%  ci(idxa,idxb) = integral(fh, 0, 1,'RelTol',0,'AbsTol',1e-6);
 end
end

 function out = eta(x)
  out = -x .* log(x) - (1-x) .* log(1-x);
 end
 function out = h(a,b)
  out = psi(0,a+b+1) - a/(a+b) * psi(0,a+1) - b/(a+b) * psi(0,b+1);
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
  subplot(2,2,1); 
  imagesc(squeeze(cf)); 
  xlabel('\alpha');ylabel('\beta')
  title(['Analytic specific heat divergence rate c(T=1)/N'])
  set(gca, 'YTick', round(linspace(1,length(as), 4)));
  set(gca, 'YTickLabel', as(round(linspace(1,length(as), 4))));
  set(gca, 'XTick', round(linspace(1,length(bs), 4)));
  set(gca, 'XTickLabel', bs(round(linspace(1,length(bs), 4))));
  set(gca, 'TickDir', 'out'); box off
  subplot(2,2,2); 
  imagesc(squeeze(cs));
  xlabel('\alpha');ylabel('\beta')
  title(['Monte Carlo results'])
  set(gca, 'YTick', round(linspace(1,length(as), 4)));
  set(gca, 'YTickLabel', as(round(linspace(1,length(as), 4))));
  set(gca, 'XTick', round(linspace(1,length(bs), 4)));
  set(gca, 'XTickLabel', bs(round(linspace(1,length(bs), 4))));
  set(gca, 'TickDir', 'out'); box off
  subplot(2,2,3); 
  imagesc(squeeze(ci)); 
  xlabel('\alpha');ylabel('\beta')
  title(['Numerical integration results'])
  set(gca, 'YTick', round(linspace(1,length(as), 4)));
  set(gca, 'YTickLabel', as(round(linspace(1,length(as), 4))));
  set(gca, 'XTick', round(linspace(1,length(bs), 4)));
  set(gca, 'XTickLabel', bs(round(linspace(1,length(bs), 4))));
  set(gca, 'TickDir', 'out'); box off
  subplot(2,2,4); 
  imagesc(squeeze(cf)./squeeze(ci)); colorbar
  xlabel('\alpha');ylabel('\beta')
  title(['Ratio of analytic and numerical integration results'])
  set(gca, 'YTick', round(linspace(1,length(as), 4)));
  set(gca, 'YTickLabel', as(round(linspace(1,length(as), 4))));
  set(gca, 'XTick', round(linspace(1,length(bs), 4)));
  set(gca, 'XTickLabel', bs(round(linspace(1,length(bs), 4))));
  set(gca, 'TickDir', 'out'); box off
%  maxRelDiff=max(vec(abs(squeeze(Ef)-squeeze(Es)))./abs(vec(squeeze(Ef))));
%  annotation('textbox', [0 0.9 1 0.1], ...
%    'String', ['Maximal relative error max\{ (E_{pred} - E_{MC}) / E_{pred} \} =  ', num2str(round(maxRelDiff*1000)/10), '%'], ...
%    'EdgeColor', 'none', ...
%    'HorizontalAlignment', 'center')  
 end
end