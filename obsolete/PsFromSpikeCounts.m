function [Ps,mynchoosek]=PsFromSpikeCounts(p_count);

N=numel(p_count)-1;

lognchoosek=gammaln(N+1)-gammaln((0:N)'+1)-gammaln(N-(0:N)'+1);
mynchoosek=round(exp(lognchoosek));

Ps=[];
for k=1:(N+1);
p=p_count(k)/(mynchoosek(k));
Ps=[Ps; repmat(p,(mynchoosek(k)),1);];
end