% Calculates the dynamic quantile test test of conditional coverage
% of Engle and Manganelli (2004).
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fDynamicQuantileTest(I,p,lags,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   p         -  Coverage rate of VaR (probability of a hit)
%   lags      -  Lags to include in test
%   sign      -  (Optional) significance level for assymptotic critical value, default 0.05
%   bootstrap -  (Optional) Indicates wheather bootstrapped p-values should
%                be returned. Takes values 'yes' or 'no', default is n
%
% OUTPUTS:
%   Test        -  Test value
%   asymptotics -  critical value of sign significance in asymptotic distribution
%   name        -  name of test
%
% Comments:     Only the lagged values of the Hitsequence is used as
%               explanatory variables. This means that the asymptotic 
%               distribution is chi(k+1) rather than chi(2k+1)
%
% EXAMPLE:
% p = 0.05;                      %Coverage rate
% T = 500;                       %Observations
% I = binornd(1,p,T,1);          %Simulates hit-sequence     
% d = fDurations(I);             %Gets durations
% fDynamicQuantileTest(I,p,10)   %Calls the "dynamic Quantile" test of Conditional coverage using 10 lag
% fDynamicQuantileTest(I,p,20)   %Calls the "dynamic Quantile" test of Conditional coverage using 20 lag
%
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fDynamicQuantileTest(I,p,lags,sign,bootstrap)

%converts hit-seq to doubles, easier for matlab mex files
if islogical(I)==1
    I=+I;
end

%checks number of input
if nargin <3
    error('Atleast 3 inputs are required.');
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
if (length(I)<2);
    error('Hit-sequence, I, is not of length >1');
end

%checks if number of lags tested is an integer
if (rem(lags,1) ~=0)
    error('Number of lags tested, k, must be a whole number.');
end

%If hits, calculate test statistic
if sum(I)~=0

    %demean hit seq
    Hit = I-p;

    %create the model matrix used in the OLS
    Z = NaN(length(I)-lags,(1+lags));
    Time = (1+lags):length(I);
    
    for i=1:lags
        if i==1;
            LagPart = Hit(Time-i);
        else 
            LagPart = [LagPart Hit(Time-i)];
        end
    end
    
    Z(Time,:)=[ones(size(I,1)-lags,1) LagPart];

    Z=Z((1+lags):length(I),1:size(Z,2));

    %Discards observations to fit with the lags
    Hit = Hit((1+lags):length(Hit),1);

    %OLS
    b = Z\Hit;

    %calculates test statistic
    Test = (b'*(Z'*Z)*b)/(p*(1-p));

else
    disp('Test set to 0, no hits in hit-sequence.');
    Test = 0;
end

%Calculats asymptotic critical value
asymptotics = chi2inv(1-sign,lags+1);

%Name of test
name = ['DQ-' num2str(lags)];

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,lags,sign,bootstrap);
end

end







