function [freq]=ps_2_count_distrib(ps,states)
%convert full probability mass function over states to count distribuion

%also see CountOnes

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
    freq(k+1)=sum(ps(theones==k));
end
end
freq=freq/sum(freq);
