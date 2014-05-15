function qIsing=qIsing(dat,H,J,Z,logit)
%probability of data dat under Ising model with H, J Z
if nargin<=3
    Z=1;
end

dat=dat*H+.5*diag(dat*J*dat');
if nargin<=4 || logit==0
qIsing=exp(dat)/Z;
else
qIsing=dat-log(Z);
end
