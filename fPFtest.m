% Calculates the Proportion of failures test statistic of Kupiec (1995) based on hit sequence
% and coverage inputs.
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fPFtest(I,p,sign,bootstrap)
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
%   asymptotics -  Critical value of sign significance in asymptotic distribution
%   name        -  Name of test
%
% EXAMPLE:
% a = 0.95;               %Coverage set to 95%
% p = 1-a;                %Coverage rate
% T = 500;                %Observations
% I = binornd(1,p,T,1);   %Simulates hit-sequence       
% fPFtest(I,p)            %Calls "Proportion of failures" test function
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fPFtest(I,p,sign,bootstrap)

%converts hit-seq to doubles, easier for matlab mex files
if islogical(I)==1
    I=+I;
end

%checks number of input
if nargin <2
    error('Atleast 2 inputs are required.');
end

%Default significance level is set to 5% for the asymptotic critical value
if exist('sign','var') == 0
    sign = 0.05;
end

%Default bootstrap to no
if exist('bootstrap','var') == 0
    bootstrap = 'no';
end

%Checks that significance level is between 0 and 1
if  (sign<=0) || (sign>=1);
    error('Significance level, p, for test is not between 0 and 1.');
end

%Checks coverage level is a of value between 0 and 1
if  (p<=0) || (p>=1);
    error('Coverage level, p, for test is not between 0 and 1.');
end

%checks hit-sequence is of length greater than 2 (a vector)
if (size(I,1)<2 || size(I,2)>1);
    error('Hit-sequence, I, is not a vector of length >1');
end

%If the observed number of violations matched the expected number, Test
%statistic is set to 0 to avoid floating point errors.
if abs(mean(I)-p)>0.000001
    %Number of observations
    n = length(I);

    %Number of failures
    x = sum(I);

    %calculates test statistic (robust version)
    if x~=0
        Test = -2*((n-x)*log(1-p)+x*log(p)-(n-x)*log(1-x/n)-x*log(x/n));
    else
        Test = -2*((n-x)*log(1-p)+x*log(p)-(n-x)*log(1-x/n));
    end
else
    Test = 0;
end

%Calculats asymptotic critical value
asymptotics = chi2inv(1-sign,1);

%Name of test
name = 'PF';

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,sign);
end

end
