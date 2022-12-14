% Calculates the discrete Duration based test of Haas (2006).
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fDurDtest(I,p,criteria,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   p         -  Coverage rate of VaR (probability of a hit)
%   criteria  -  Choosing between testing conditional coverage (cc) 
%                or independence (ind). Takes values 'cc' or 'ind'
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
% a = 0.95;                     %Coverage set to 955%
% p =1-a;                       %Coverage rate
% T = 500;                      %Observations
% I = binornd(1,p,T,1);         %Simulates hit-sequence     
% fDurDtest(I,p,'cc')           %calls the "discrete duration" test of cc
% fDurDtest(I,p,'ind')          %calls the "discrete duration" test of ind
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fDurDtest(I,p,criteria,sign,bootstrap)

%converts hit-seq to doubles, easier for matlab mex files
if islogical(I)==1
    I=+I;
end

%checks number of input
if nargin <3
    error('Atleast 3 inputs are required.');
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

% Degrees of freedom for asymptotic distribution is set
if strcmp(criteria,'cc')
    df = 2;
elseif strcmp(criteria,'ind')
   df = 1;
else
    error('Incorrect input, you must specify either cc or ind criteria');
end

% create censoring series
c = zeros(size(D,2),1);
if I(1)==0; c(1) = 1; end
if I(size(I,1))==0; c(size(D,2)) = 1; end

%checks if enough durations are used, else test set to 0
if sum(I)>2 || (sum(I)==2 && (I(1)+I(end))<2)
    
    % Discrete weibull is estimated and logL found
    parms = wblDiscretefit(D,c);
    loglW = fDscreteWeibulllogl(D,c,parms);

    % Geometric is estimated/fixed based on what criteria is tested and logL is found
    % degrees of freedom for asymptotic distribution is set
    if strcmp(criteria,'cc')
        a = -log(1-p);
    elseif strcmp(criteria,'ind')
       a = wblDiscretefit2(D,c); 
    end
    parms = [a 1];
    loglG = fDscreteWeibulllogl(D,c,parms);

    % calculates test statistic
    Test = -2*loglG+2*loglW;

elseif sum(I)<2
    disp('Not enough durations, test set to 0.');
    Test =0;
end

% Calculats asymptotic critical value
asymptotics = chi2inv(1-sign,df);

% Name of test
name = ['DurDtest-' criteria];

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,criteria,sign,bootstrap);
end

end

%%
% Estimates the discrete Weibull distribution of Najagawa and Osaki (1975)
% in its transformed version (more robust)
%
% USAGE:
%   Estimates = wblDiscretefit(D,c)
%
% INPUTS:
%   D        -  Durations
%   c        -  censoring series
%
% OUTPUTS:
%   Estimates -  Estimates
%
% Author:       Thor Nielsen
% Date:         14/08/2013

function Estimates = wblDiscretefit(D,c)

%Set nummerical stuff
InitialGuess = [0.10 1];
LowerBound = [0 0];

%Maximum likelihood estimates are found
options = optimoptions('fmincon','Algorithm','sqp','Display','off','Diagnostics','off');
[Estimates,~,~,~,~,~,~] = fmincon(@(parms) -fDscreteWeibulllogl(D,c,parms),InitialGuess,[],[],[],[],LowerBound,[],[],options);

end

%%
% Logl for the discrete Weibull distribution of Najagawa and Osaki (1975)
% in its transformed version (more robust)
%
% USAGE:
%   loglW = fDscreteWeibulllogl(D,c,parms)
%
% INPUTS:
%   D        -  Durations
%   c        -  censoring series
%   parms    -  Vector of parameters
%
% OUTPUTS:
%   loglW     - logL
%
% Author:       Thor Nielsen
% Date:         14/08/2013

function loglW = fDscreteWeibulllogl(D,c,parms)
%Gets parameters
a = parms(1);
b = parms(2);

% Calculates logl
loglW = c(1)*log(exp(-a^b*(D(1)-1)^b))+(1-c(1))*log(exp(-a^b.*(D(1)-1)^b)-exp(-a^b*D(1)^b));
loglW = loglW+c(end)*log(exp(-a^b*(D(end)-1)^b))+(1-c(end))*log(exp(-a^b*(D(end)-1)^b)-exp(-a^b*D(end)^b));
for i=2:(length(D)-1)
 loglW=loglW+log(exp(-a.^b.*(D(i)-1).^b)-exp(-a.^b*D(i).^b));
end

end

%%
% Estimates the geometric distribution 
%
% USAGE:
%   Estimates = wblDiscretefit2(D,c)
%
% INPUTS:
%   D        -  Durations
%   c        -  censoring series
%
% OUTPUTS:
%   Estimates -  Estimates
%
% Author:       Thor Nielsen
% Date:         14/08/2013

function Estimates = wblDiscretefit2(D,c)

%Set nummerical stuff
InitialGuess = 0.10;
LowerBound = 0;

%Calls nummerical stuff
options = optimoptions('fmincon','Algorithm','sqp','Display','off','Diagnostics','off');
[Estimates,~,~,~,~,~,~] = fmincon(@(a) -fDscreteWeibulllogl(D,c,[a 1]),InitialGuess,[],[],[],[],LowerBound,[],[],options);

end
