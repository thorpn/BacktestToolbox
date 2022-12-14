% Calculates the Markov test of Christoffersen (1998) based on hit sequence
% and coverage inputs.
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fMarkovtest(I,p,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   p         -  Coverage rate of VaR (probability of a hit)
%   sign      -  (Optional) significance level for assymptotic critical value, default 0.05
%   bootstrap -  (Optional) Indicates wheather bootstrapped p-values should
%                be returned. Takes values 'yes' or 'no', default is n
%
% OUTPUTS:
%   Test        -  Test value
%   asymptotics -  critical value of sign significance in asymptotic distribution
%   name        -  name of test
%
% EXAMPLE:
% a = 0.95;                 %Coverage set to 955%
% p = 1-a;                  %Coverage rate
% T = 500;                  %Observations
% I = binornd(1,p,T,1);     %Simulates hit-sequence       
% fMarkovtest(I,p)          %Calls "Markov" test function
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% Date:     15-04-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fMarkovtest(I,p,sign,bootstrap)

%converts hit-seq to doubles, easier for matlab mex files
if islogical(I)==1
    I=+I;
end

%checks number of input
if nargin <2
    error('Atleast 2 inputs are required.');
end

%Default bootstrap to no
if exist('bootstrap','var') == 0
    bootstrap = 'no';
end

% %checks number of observations and coverage for nummerical issues
% if (length(I)>3000 && p<0.010001) || length(I)>10000 ;
%     error('Too many observations, calculate the test statistic using the log rules shown in the paper for nummerical robustness.');
% end

%Default significance level is set to 5% for the asymptotic critical value
if exist('sign','var') == 0
    sign = 0.05;
end

%Checks that significance level is between 0 and 1
if  (sign<=0) || (sign>=1);
    error('Significance level, p, for test is not between 0 and 1.');
end

%checks hit-sequence is of length greater than 2 (a vector)
if (length(I)<2);
    error('Hit-sequence, I, is not of sufficient length');
end

%Calls the generalized test, with state = 1 this is the Christoffersen test
[Test, asymptotics, name] = fGeneralizedMarkovtest(I,p,1,'ind',sign,bootstrap);

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,1,'ind',sign,bootstrap);
end

end


