function checkProperTemperatureRange(i, nT, Tmin, Tmax, nRep, nSamples, burnIn)
% this function serves to test some temperature magic on the lambdaTrue
% used for the previous test runs for the iterative scaling code

if nargin < 7
    burnIn = 32000;
end
if nargin < 6
    nSamples = 32000;
end
if nargin < 5
    nRep = 30;
end
if nargin < 4
    Tmax = 2;
end
if nargin < 3
    Tmin = 0.5;
end
if nargin < 2
    nT = 21;
end

addpath(genpath('/home/marcel/criticalityIterScaling/code/'))% get codebase

load('/home/marcel/criticalityIterScaling/s1_res_small.mat')
lambdaTrue = out.lambdaTrue; clear out; % get so far used lambdaTrue

Ts = linspace(Tmin, Tmax, nT); % this is the same for all calls with all i

x0 = 100; % first run will have to initialized SOMEHOW 

tempTracesEfy = zeros(5151,nRep); % output variable

lambdaT = lambdaTrue / (Ts(i)); % parameter vector for this temperature run
 
for j = 1:nRep
 disp(['checking draw #', num2str(j),' out of ', num2str(nRep)])
 [tempTracesEfy(:,j),~,x0] = maxEnt_gibbs_pair_C(nSamples, burnIn, ...
       lambdaT, x0, 'cluster');
end
save(['/home/marcel/criticalityIterScaling/tempTest', num2str(i),'of',...
      num2str(nT),'.mat'], 'tempTracesEfy', 'lambdaT', 'Ts');
end
