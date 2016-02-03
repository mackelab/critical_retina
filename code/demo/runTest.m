function runTest(fname, nSample, a, tau, maxIter, ifVK, ifbwVK)

if nargin < 6
    ifVK = true;
end

if nargin < 7
    ifbwVK = true;
end

ifSave = true;

load('/home/marcel/criticalityIterScaling/data/testDataBest.mat')
testPackageCluster(fname,100,32000,32000,1000, nSample, 200, maxIter, 1, ...
                               [], [], lambdaTrue ,mfxTrain, ...
                               a, tau, ifSave, ifVK, ifbwVK)
                           
end