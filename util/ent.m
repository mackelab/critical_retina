function ent=ent(P,logbase)
%calculate entropy given probability vector, for basis e or 2. If no
%argument is given, uses base E.

P=P/sum(P);

P(P<exp(-700))=exp(-700);

if nargin==1
    logP=log(P);
elseif logbase==2
    logP=log2(P);
end


logP(isnan(logP))=0;

ent=-sum(P.*logP);
