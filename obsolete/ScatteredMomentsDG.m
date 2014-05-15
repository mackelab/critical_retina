%run simulations for within-model class bias of Ising models:

if 1
    clear all
    addpath([pwd,'/../../FitMaxEnt']);
    addpath ~/CopyMPI/Code/minFunc/;
    addpath ~/CopyMPI/Code/minFunc/;
    addpath ~/CopyMPI/MaxEnt/matlab/MaxEntBias/functions
    addpath ~/CopyMPI/BinCorr/matlab/FlatModels
    addpath ~/CopyMPI/BinCorr/matlab/FlatDG
    addpath ~/CopyMPI/BinCorr/matlab/functions/
    addpath ~/CopyMPI/BinCorr/matlab/DGfitting/
    basename='~/CopyMPI/MaxEnt/Data/ScatteredMomentsDG';
    thisdir=pwd;
    cd('~/CopyMPI/matlab/functions'), pathchoices=4; startup, cd(thisdir)
    
    mkdir(basename)
    callfile=mfilename;
    dimo=[5]; %use Ising model with 2, 3 5, 10, 15 dimensions
    
    %means=[.05,.1,.5]; %try out means 0.01, .1, .5
    means=norminv([.05:.05:.5]);
    startcorrs=[0:.05:.5];
    scatterstds=[0:.05:.5]; %standard devation of Js for population size 1, for population size n, std is divided by sqrt(n)
    numrethrowparams=25; %for each experiment, try 5 random throws of parameters J
    Ns=[10:10:200,round(10.^linspace(2,5,25))]; %number of 'datapoints' in simulation.
    Ns=unique(sort(Ns));
    
    numrepeats=10000;
    
    
    
    
    %for each setting of parameters, create many data-sets to average over, so that we can calculate the expected bias, and the variance of the estimator
    numrepeats_save=5; %save detailed results for the first 10 data-sets:
    minmean=1e-20; %minimum allowed mean, any smaller mean is trunated at this value (otherwise maxent model is not defined, or lives in lower-dimensional space)
    
    %options for parameter learning:
    fitoptions.TolX=1e-50;
    fitoptions.TolFun=1e-50;
    fitoptions.display='off';
    fitoptions.MaxIter=5000;
    fitoptions.restarts=2;
    [features,description,x]=SetupFeaturesMaxEnt(dimo,2);
    
    %first, set up all the parameters and save them:
    savenameparams=[basename,'Params'];
    
    for k=1:numel(means);
        k
        for kk=1:numel(startcorrs)
            kk
            for kkk=1:numel(scatterstds);
                tic
                for kkkk=1:numrethrowparams
                    
                    p.fmeanorig=means(k);
                    p.startcorr=startcorrs(kk);
                    p.dim=dimo;
                    p.scatterstd=scatterstds(kkk);
                    p.rethrow=kkkk;
                    
                    %first, make random gamma vector and correlation matrix:
                    
                    Lambda0=(eye(p.dim)*(1-p.startcorr)+ones(p.dim)*p.startcorr);
                    G0=chol(Lambda0);
                    G=G0+triu(randn(size(G0)),1)*p.scatterstd;
                    p.Lambda=Cov2Corr(G'*G);
                    gamma0=ones(p.dim,1)*p.fmeanorig;
                    
                    p.gamma=gamma0+randn(p.dim,1)*p.scatterstd;
                    [p.mu, p.Rho]=Lambda2Rho(p.gamma, p.Lambda);
                    p.meanmean=mean(p.mu);
                    p.medmean=median(p.mu);
                    p.harmman=exp(mean(log(p.mu)));
                    p.stdmean=std(p.mu);
                    p.rangemean=max(p.mu)-min(p.mu);
                    
                    corrs=p.Rho(triu(ones(size(p.Rho)),1)==1);
                    p.meancorr=mean(corrs);
                    p.medcorr=median(corrs);
                    p.harmcorr=exp(mean(log(corrs)));
                    p.stdcorr=std(corrs);
                    p.rangecorr=max(corrs)-min(corrs);
                    
                    
                    %DG first:
                    p.P_DG=ProbsDG(x,p.gamma,p.Lambda)';
                    mean_features_DG=features'*p.P_DG;
                    p.Entropy_DG=Ent(p.P_DG);
                    
                    %then Ising model
                    [p.lambda_q2,P.logZ_q2, P_q2, p.mean_features_q2,output]=FitMaxEntLinear(features,mean_features_DG', fitoptions);
                    P_q2=exp(P_q2);
                    p.Entropy_q2=ent(P_q2);
                    
                    
                    %then Bias calculations
                    Sigma_q2=wcov(features,P_q2);
                    p.Bx=CalcBx(features, p.mean_features_q2,Sigma_q2);
                    p.Bias_q2=p.Bx'*P_q2;
                    p.Bias_DG=p.Bx'*p.P_DG;
                    featuresBiased=[features,p.Bx];
                    [p.bDash,Hx,Bx,VarB,EBdeltag,SigmaBdeltag]=CalcbDashAtZero(features,p.mean_features_q2, Sigma_q2, P_q2);
                    
                    %then MaxBias model
                    meansBiased=[p.mean_features_q2,p.Bias_DG];
                    try
                        [p.lambda_biased,p.logZ_biased, P_biased, p.means_biased,output]=FitMaxEntLinear(featuresBiased,meansBiased, fitoptions);
                        P_biased=exp(P_biased);
                        p.Accuracy=max(abs(p.means_biased-meansBiased));
                        p.Entropy_biased=ent(P_biased);
                    catch
                       % keyboard
                        p.lambda_biased=nan;
                        p.logZ_biased=nan;
                        P_biased=nan;
                        p.means_biased=nan;
                        p.P_biased=nan;
                        %p.DeltaS_DG=nan;
                        %P.DeltaS_min=nan;
                        %p.DeltaS_perturb=nan;
                        p.Accuracy=nan;
                    end
                    
                    
                   % p.started=nan;
                   % p.finished=nan;
                    %p.filename=sprintf('%s/Results_%03d_%03d_%03d_%03d',basename,k,kk,kkk,kkkk);
                    params(k,kk,kkk,kkkk)=p;
                    %need checker that tells us when there are problems!!!!!
                    
                end
                toc
                save([basename,'/HQ']);
            end
        end
    end
else
    %clear all
    thisdir=pwd;
    cd('~/CopyMPI/matlab/functions'), pathchoices=4; startup, cd(thisdir)
    addpath([pwd,'/../../FitMaxEnt']);
    addpath ~/CopyMPI/Code/minFunc/;
    basename='~/CopyMPI/MaxEnt/Data/BiasIsingScattered';
    load([basename,'/HQ']);
end
