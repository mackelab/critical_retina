function [lambda,logZ,logP,fitmeans,output]=fit_maxent_mpf(x, fitoptions)

[n, d] = size( x ); % watch out, n and d will switch meaning when packing
                    % data x into data.x for the MPF code. Sorry for that.

%if starting point is provided, use it:
if isfield(fitoptions,'lambda0')
    lambda=fitoptions.lambda0(:);
else
    %otherwise use zeros:
    lambda=zeros( d + d*(d-1)/2 + (d+1), 1);
               %  h +     J     +   L 
end

 % Parameter convention-conversion fun.
  J = diag(lambda(1:d));                  % fiddle contents of h into J
  J(logical(triu(ones(size(J)),1))) = ... % fiddle contents of J into J
                              lambda(d+(1:d*(d-1)/2)); 
  J = (J + J')/2; % MPF actually requires some real matrix computations of
                  % J with the data matrix x, thus it appears a bit tedious
                  % to work with upper-triangular J and pack/unpack each 
                  % call and to instead simply work with n^2 entries for J. 
    
  L = lambda(end-d:end); % assumes d+1 entries, i.e. a feature for K = 0

  lambda = [J(:); L(:)]; % take note, this is what the following MPF code 
                         % assumes the parameters to be structured as
   
 % Data preprocessing (data does not change over minFunc calls, do it once)
  data.x      = x';
  data.counts = sum(data.x,1);
  data.mask   = ones(d,n);    
  for i = 1:n
    idx = find(sum(data.x ~= data.x(:,i)*ones(1,n))==1);
    for j = 1:length(idx)
     data.mask(data.x(:,i)~=data.x(:,idx(j)),i) = 0; 
    end
  end
 % Numerical optimization step
  [lambda,f,exitflag,output] = minFunc( @K_dK_ising_PK, lambda, fitoptions, data );
  output.fs=f;
  output.exitflag=exitflag;
  
 % More parameter convention-conversion fun.  
  J = reshape(lambda(1:d^2),d,d);
  h = diag(J);
  J = 2*J(logical(triu(ones(size(J)),1)));
  L = lambda(end-d+1:end);                     % CURRENTLY JUST DROPPING THE PARAMETER FOR K = 0 !
  lambda = [ h(:); J(:); L(:) ];
  
   weights = 0; % currently cannot set them otherwise with this method. 
   
   fx = setup_features_maxent(x, 'ising_count');
  [logP,logZ,~,fitmeans]= logPMaxEnt(fx,lambda,[],weights);
  
end
  
  