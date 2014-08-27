
function [xSampled,p1s,p0s,ks] = maxEnt_gibbs(nSamples, burnIn, lambda, x0, model)

p1s = zeros(burnIn+nSamples+1,1);
p0s = zeros(burnIn+nSamples+1,1);
ks  = zeros(burnIn+nSamples+1,1);

d = length(x0);

xSampled = zeros(d, burnIn+nSamples+1); % +1 for x0
xSampled(:,1) = x0;
fxc = maxEnt_features(x0, model); 

m   = zeros(length(fxc),d); % terms corresponding to features for h
for k = 1:d % compute masks
   m(k,k)= 1;                  %   
   tmp = zeros(d,d);   % add terms corresponding to features for J
   tmp(k,:) = 1; tmp(:,k) = 1;
   tmp = tmp(logical(triu(ones(d,d),1)));
   m(d+(1:d*(d-1)/2),k)   = tmp(:);        
end

for i = 2:nSamples+burnIn        
  k = randi(d);
        
  xc = xSampled(:,i-1); % current sample       
        
 
  % Compute probability
   v = m(:,k); % current mask for choosing correct feature components. 
               % Last entries depend on current sample xc for 
               % model = 'ising_count', 'ising_count_l_0'

   idl = sum(xc) - xc(k);
   switch model
       case 'ising'
        p0 = 1;              
       case 'ising_count'
        v(end-d+1:end) = fxc(end-d+1:end); % last d entries are indicators
        if idl>0, p0 = exp(-lambda(end-d+1+idl)); else p0 = 1; end
       case 'ising_count_l_0'
        v(end-d:end) = fxc(end-d:end);   % last d+1 entries are indicators
        p0 = exp(-lambda(end-d+idl));
   end
   xc1 = xc; xc1(k) = 1;
   fxc1 = maxEnt_features(xc1, model); % features, needed for v for V(K)
   p1 = exp(-(lambda.*v)' * fxc1);

   p1s(i) = p1;       
   p0s(i) = p0;
   ks(i) = k;

   p1  = p1 / (p0 + p1); % normalization step
 
  % Draw new k-th entry
    xk = rand(1) < p1;
  
  % Update chain  
    xSampled(:,i) = xSampled(:,i-1);
    xSampled(k,i) = xk;
end
   % Discard burnin samples
 xSampled = xSampled(:,burnIn+2:end);
 
end

function fx = maxEnt_features(x, model)
 d = length(x);
 switch model
    case {2,'ising'}
        pairs=nchoosek(1:d,2);
        %description=[[1:d; (1:d)*nan],pairs];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2))];
    case 3
        pairs=nchoosek(1:d,2);
        triplets=nchoosek(1:d,3);
        %description=[[1:d; (1:d)*nan;(1:d)*nan],[pairs; pairs(1,:)*nan],triplets];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));x(pairs(:,1)).*x(pairs(:,2)).*x(triplets(:,3))];
    case 'ising_count'
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d,1);
        sum_x = sum(x);
        if sum_x>1
         count_indicators(sum_x)=1;
        end
        %description=[[1:d; (1:d)*nan],pairs,[1:d; (1:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
        %keyboard
    case 'ising_count_l_0' % same as 'ising_count', but has feature for K=0
        pairs=nchoosek(1:d,2);
        count_indicators=zeros(d+1,1);
        count_indicators(sum(x)+1) = 1; 
        %description=[[1:d; (1:d)*nan],pairs,[0:d; (0:d)*nan]];  
        fx=[x;x(pairs(:,1)).*x(pairs(:,2));count_indicators];
end

end
