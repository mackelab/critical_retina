function [count_distribs,Zs,models]=flat_model_trace_temp(count_distrib,betas);
%fiven a model with specificed spike-count distribution, calculate
%spike-count distributions for different inverse temperatures beta;
%
%formula: for states x, P_b(x)= 1/Z_b exp(beta*lambda'*F(x));
%for spike-counts k, P_b(k)= 1/Z_b  (N choose k)^{1-beta}   P_1(k);



for k=1:numel(betas);
    [count_distribs(k,:),Zs(k,:),models(k)]=flat_model_vary_temp(count_distrib,betas(k));
    model(k).beta=betas(k);
end
