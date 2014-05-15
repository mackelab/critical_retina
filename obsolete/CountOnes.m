function [freq]=CountOnes(states,P)
%count the number of ones in patterns supplied in "states";
if ~islogical(states)
    states=DecToBinary(states);
end

theones=sum(states,2);

if nargin==1
    for k=0:size(states,2)
        freq(k+1)=sum(theones==k);
    end
    
else
for k=0:size(states,2)
    freq(k+1)=sum(P(theones==k));
end
end
freq=freq/sum(freq);
