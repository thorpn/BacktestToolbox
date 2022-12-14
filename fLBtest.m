% Calculates the Ljung-Box test statistic based on hit sequence
% and coverage inputs.
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fLBtest(I,p,k,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   k         -  Number of lags tested
%   sign      -  (Optional) significance level for assymptotic critical value, default 0.05
%   p         -  Coverage rate of VaR (probability of a hit)
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
% fLBtest(I,p,1)            %Calls "Ljung-Box" test function examining 1 correlations
% fLBtest(I,p,10)           %Calls "Ljung-Box" test function examining 10 correlations
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fLBtest(I,p,k,sign,bootstrap)

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

%Default significance level is set to 5% for the asymptotic critical value
if exist('sign','var') == 0
    sign = 0.05;
end

%Checks that significance level is between 0 and 1
if  (sign<=0) || (sign>=1);
    error('Significance level, p, for test is not between 0 and 1.');
end

%checks if number of autcorrelations tested is an integer
if (rem(k,1) ~=0)
    error('Number of autcorrelations tested, k, must be a whole number.');
end

%checks hit-sequence is of length greater than 2 (a vector)
if (length(I)<2);
    error('Hit-sequence, I, is not of length >1');
end

if sum(I)~=0;
  
    T = length(I);

    %Gets autocorrelations
    ACF = autocorr(I,k);

    %Removes lag 0
    ACF = ACF(2:end); 

    %calculates test statistic
    stat = T*(T+2)*cumsum((ACF.^2)./(T-(1:k))');
    Test = stat(k);


    %Calculats asymptotic critical value
    asymptotics = chi2inv(1-sign,k);

    %Name of test
    name = ['LB-' num2str(k)];

    %Calculates the bootstrapped p-value
    if strcmp(bootstrap,'yes')
        varargout{1} = fBootPval(name,Test,I,p,k,sign,bootstrap);
    end
    
else
    
    Test = 0;

end

%Calculats asymptotic critical value
asymptotics = chi2inv(1-sign,k);

%Name of test
name = ['LB-' num2str(k)];

end






