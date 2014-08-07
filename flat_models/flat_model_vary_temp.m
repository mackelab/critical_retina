function [count_distrib,Z,model]=flat_model_vary_temp(count_distrib,beta);
%fiven a model with specificed spike-count distribution, calculate
%spike-count distributions for different inverse temperatures beta;
%
%formula: for states x, P_b(x)= 1/Z_b exp(beta*lambda'*F(x));
%for spike-counts k, P_b(k)= 1/Z_b  (N choose k)^{1-beta}   P_1(k);

N=numel(count_distrib)-1;
lognchoosek=gammaln(N+1)-gammaln((0:N)+1)-gammaln(N-(0:N)+1);

logP=log(count_distrib);

count_distrib=zeros(1,numel(count_distrib));


logP_beta= beta*logP+ (1-beta)* lognchoosek;
P_beta=exp(logP_beta);
Z=sum(P_beta);
count_distrib=P_beta/Z;

if nargout==3
    model.N=N;
[model.meancount,model.varcount]=calc_mean_var(count_distrib,0:N);
[model.mean,model.corr]=meanvar_count_2_meancorr(model.meancount,model.varcount,N);
[model.entropy,model.entropy_count]=entropy_flat_model(count_distrib);
[model.var_log_probs,model.mean_log_probs]= flat_model_var_log_probs(count_distrib,exp(1));

model.count_distrib=count_distrib;
%model.patterns=exp(lognchoosek);
end
