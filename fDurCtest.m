% Calculates the Duration based test of independence (Christoffersen et al. 2004)
%
% USAGE:
%   [Test, asymptotics, name] = fDurCtest(I,p,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   D         -  Duration sequence, D, in a column vector
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
% a = 0.95;                     %Coverage set to 955%
% p = 1-a;                      %Coverage rate
% T = 500;                      %Observations
% I = binornd(1,p,T,1);         %Simulates hit-sequence
% fDurCtest(I,p)     %calls the "continuous duration" test of cc
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name] = fDurCtest(I,p,sign,bootstrap)

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

% Checks if toolbox is available, then gets durations
if exist('fDurations') == 2
    %Gets Durations and censored
    D = cell2mat(fDurations(I));
else
    error('fDurations not in memory, add the BackTestToolbox before using this test function.');
end

% Default significance level is set to 5% for the asymptotic critical value
if exist('sign','var') == 0
    sign = 0.05;
end

% Checks that significance level is between 0 and 1
if  (sign<=0) || (sign>=1);
    error('Significance level, p, for test is not between 0 and 1.');
end

% checks duration-sequence is a double
if (~isa(D,'double'));
    error('Durations, d, is not a double.');
end

% checks hit-sequence is of length greater than 2 (a vector)
if (length(I)<2);
    error('Hit-sequence, I, is not of length >1');
end

% create censoring series
c = zeros(size(D,2),1);
if I(1)==0; c(1) = 1; end
if I(size(I,1))==0; c(size(D,2)) = 1; end

%checks if enough durations are used, else test set to 0
if sum(I)>2 || (sum(I)==2 && (I(1)+I(end))<2)

    %%%%%%% Weibull is fitted and logL found
    parmhat = wblfit(D);                                                %Fits the Weibull distribution
    f = @(d) wblpdf(d,parmhat(1),parmhat(2));                           %pdf
    S =  @(d) 1-wblcdf(d,parmhat(1),parmhat(2));                        %survivor function

    loglW = c(1)*log(S(D(1)))+(1-c(1))*log(f(D(1)));                    %first observation LogL
    loglW = loglW+c(end)*log(S(D(end)))+(1-c(end))*log(f(D(end)));      %last observation LogL
    loglW = loglW+sum(log(f(D(2:(length(D)-1)))));                      %rest of LogL

    %%%%%%% Exponential is fitted and logL found
    parmhat = expfit(D);                                                %Fits the Exponential distribution
    f = @(d) exppdf(d,parmhat);                                         %pdf
    S =  @(d) 1-expcdf(d,parmhat);                                      %survivor function

    loglE = c(1)*log(S(D(1)))+(1-c(1))*log(f(D(1)));                    %first observation LogL
    loglE = loglE+c(end)*log(S(D(end)))+(1-c(end))*log(f(D(end)));      %last observation LogL
    loglE = loglE+sum(log(f(D(2:(length(D)-1)))));                      %rest of LogL

    %calculates test statistic
    Test = -2*loglE+2*loglW;

elseif sum(I)<2
    disp('Not enough durations, test set to 0.');
    Test =0;
end
    
%Calculats asymptotic critical value
asymptotics = NaN;

%Name of test
name = 'DurCtest';

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,sign,bootstrap);
end

end

