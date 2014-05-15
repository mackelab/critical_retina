function [gamma,lambda,count_distrib,model]=fit_flat_dg(mu,rho,N);
%function [gamma,lambda]=fit_flag_dag(mucount,varcount);
%
%fit flat dg model to population data with given mean count and variance of
%counts. 
%
%
%
%
[mucount,varcount]=meancorr_2_meanvar_count(mu,rho,N);
[gamma,lambda]=Rho_2_Lambda(mu,rho);


count_distrib=flat_dg_count_distrib(gamma,lambda,[0:N],N);

if nargout==4
    model.mu=mu;
    model.rho=rho;
    
    model.gamma=gamma;
    model.lambda=lambda;
    model.mucount=mucount;
    model.varcount=varcount;
    model.N=N;
    model.count_distrib=count_distrib;
    [model.entropy,model.entropy_count]=entropy_flat_model(count_distrib,2);
end
