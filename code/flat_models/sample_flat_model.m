function [s,hist_s]=sample_flat_model(pcount,nsamples)
%sample from a 'flat' model given the distribution
%of spike counts. 

n=numel(pcount)-1;
pcount=pcount/sum(pcount);
cum_p=[0,cumsum(pcount)];

s=false(n,nsamples);
ps=rand(nsamples,1);
%try
[hist_s,count_picks]=histc(ps,cum_p);
%catch
%    keyboard
%end
hist_s=hist_s(1:end-1);
count_picks=count_picks-1;
for k=1:nsamples
    s(randsample(n,count_picks(k)),k)=true;
  %  keyboard
end
%s=s==1;
%keyboard
