function runTest(fname, nSample, a, tau, maxIter)

load('/home/marcel/criticalityIterScaling/data/testDataBestHard.mat')
testPackageCluster(fname,100,32000,32000,1000, nSample, 200, maxIter, 1, ...
                               [], [], lambdaTrue ,mfxTrain, ...
                               a, tau)
                           
end