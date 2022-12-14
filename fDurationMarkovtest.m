% Calculates the Markov duration test of Pajhede (2014).
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fDurationMarkovtest(I,p,lags,criteria,sign,bootstrap)
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
% a = 0.95;                                        %Coverage set to 955%
% p = 1-a;                                         %Coverage rate
% T = 500;                                         %Observations
% I = binornd(1,p,T,1);                            %Simulates hit-sequence       
% fDurM(I,p,20,'cc')                               %calls Markov duration test of cc using 20 lags 
% fDurM(I,p,20,'ind')                              %calls Markov duration test of ind using 20 lags 
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fDurationMarkovtest(I,p,lags,criteria,sign,bootstrap)

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
   df = lags+1;
elseif strcmp(criteria,'ind')
   df = lags;
else
    error('Incorrect input, you must specify either cc or ind criteria');
end

%Allocates memory for counts
P = zeros(lags,1);

if sum(I(1:end))~=0

    %Counts n's
%     [N11, N10, N00, N01] = fCountHitsDuration(I,states);
    [N11, N10, N00, N01] = fCountHitsDuration_mex(I,lags);
    
    %Estimates transition probabilities
    P=N11./(N11+N10);
    pi01 = N01/(N00+N01);
    
    %Alternative p is based on what criteria is tested
    if strcmp(criteria,'cc')
       pi2 = p; 
    elseif strcmp(criteria,'ind')
       pi2 = mean(I(lags+1:end));
    end
    
    %%calculates test statistic (robust)    
    Parts = -log(1-P(:)).*N10(:)-log(P(:)).*N11(:);
    Test = -2*(log(1-pi2)*(N00+sum(N10))+nansum(log(pi2)*(sum(N11)+N01))-nansum(log(1-pi01)*N00)-nansum(log(pi01)*N01)+nansum(Parts));
else
    disp('Test set to 0, no hits in hit-sequence.');
    Test = 0;
end

%Calculats asymptotic critical value
asymptotics = chi2inv(1-sign,df);

%Name of test
name = ['DurM-' criteria '-' num2str(lags)];

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,lags,criteria,sign,bootstrap);
end

end