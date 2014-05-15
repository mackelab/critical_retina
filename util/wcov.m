function [c,m]=wcov(x,P);
%wheighted covariance matrix and mean, i.e. mean and covariance of x under
%the distribution P

wx=bsxfun(@times,full(x),P);
m=sum(wx,1);

c=x'*wx-m'*m;
