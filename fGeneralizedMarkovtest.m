% Calculates a more general version of the Markov test of Pajhede (2014)
%
% USAGE
%   [Test, asymptotics, name, varargout] = fGeneralizedMarkovtest(I,p,lags,criteria,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   lags    -  Number of lags against which we want power
%   p         -  Coverage rate of VaR (probability of a hit)
%   criteria  -  Choosing between testing conditional coverage (cc)
%                or independence (ind). Takes values 'cc', 'uc' or 'ind'
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
% a = 0.95;                               %Coverage set to 95%
% p = 1-a;                                %Coverage rate
% T = 500;                                %Observations
% I = binornd(1,p,T,1);                   %Simulates hit-sequence       
% fGeneralizedMarkovtest(I,p,20,'cc')     %Calls "generalized Markov" test of CC function using 20 lags 
% fGeneralizedMarkovtest(I,p,20,'ind')  %Calls "generalized Markov" test of Ind function using 20 lags 
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     14-08-2014
% Version:  1.0
%%

function [Test, asymptotics, name, varargout] = fGeneralizedMarkovtest(I,p,lags,criteria,sign,bootstrap)

%converts hit-seq to doubles, easier for matlab mex files
if islogical(I)==1
    I=+I;
end

%checks number of input
if nargin <4
    error('Atleast 4 inputs are required.');
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

%checks hit-sequence is of length greater than 2 (a vector)
if (length(I)<(lags+1));
    error('Hit-sequence, I, is not of sufficient length');
end

%Degrees of freedom for asymptotic distribution
if strcmp(criteria,'cc')
   df = 2;
elseif strcmp(criteria,'ind')
   df = 1;
else
    error('Incorrect input, you must specify either cc or ind criteria');
end
   
if sum(I(1:end))~=0

    %Count n's
%   [n00, n01, n10, n11]  = fCountHitsGeneralized(I,states);
    [n00, n01, n10, n11]  = fCountHitsGeneralized_mex(I,lags);
    
    %Transition probability estimates are found
    pi01 = (n01)/(n00+n01);
    pi11 = (n11)/(n10+n11);   
    
    %Alternative p is based on what criteria is tested
    if strcmp(criteria,'cc')
       pi2 = p; 
    elseif strcmp(criteria,'ind')
       pi2 = (n01+n11)/(n00+n10+n01+n11);
    end
    
%     %calculates test statistic (non-robust)
%     L0 = (1-pi2)^(n00+n10)*pi2^(n01+n11);
%     L1 = (1-pi01)^n00*pi01^n01*(1-pi11)^n10*pi11^n11;
%     Test = -2*(log(L0)-log(L1));

    %calculates test statistic (robust)
    if pi11 == 0
        Test = -2*(log(1-pi2)*(n00+n10)+nansum(log(pi2)*(n01+n11)-log(1-pi01)*n00-log(pi01)*n01));
    else
        Test = -2*(log(1-pi2)*(n00+n10)+nansum(log(pi2)*(n01+n11)-log(1-pi01)*n00-log(pi01)*n01-log(1-pi11)*n10-log(pi11)*n11));
    end
    
else
%     disp('Test set to 0, no hits in hit-sequence.');
    Test = 0;
        
end

%Calculats asymptotic critical value
asymptotics = chi2inv(1-sign,df);

%Name of test
name = ['Markov-' criteria '-' num2str(lags)];

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,lags,criteria,sign,bootstrap);
end

end



