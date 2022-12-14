

%Add toolbox to workspace (insert the path to where you have the toolbox)
addpath(genpath('insert path to toolbox here'));

%Simulate a hit-sequence
p = 0.05;                               %Coverage rate calculated to 5% value-at-risk
T = 500;                                %Number of observations
I = binornd(1,p,T,1);                   %Simulates hit-sequence       

%Calls "generalized Markov" test of CC function using 10 lags 
Test = fGeneralizedMarkovtest(I,p,10,'cc');

%Uses bootstrap of Dufour to get p-value at and sets size to 5%
[Test, asymptotics, name, pval] = fGeneralizedMarkovtest(I,p,10,'cc',0.05,'yes'); 

%Presents results
disp(['Test Name:  ' name]);
disp(['Test statistic = ' num2str(Test)]);
disp(['Asymptotic 95% Critical value = ' num2str(asymptotics)]);
disp(['Asymptotic p-value  = ' num2str(1-chi2cdf(Test,2))]);
disp(['Bootstrapped p-value  = ' num2str(pval)]);


