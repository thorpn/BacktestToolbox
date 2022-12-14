% Calculates the GMM test of Candelon et. al. (2008)
%
% USAGE:
%   [Test, asymptotics, name, varargout] = fGMMtest(I,p,k,criteria,sign,bootstrap)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   p         -  Coverage rate of VaR (probability of a hit)
%   k         -  Moments used
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
% Comments:  When k=1 only the unconditional coverage is tested, authors
%            Authors found optimal power using k=3 when testing conditional
%            coverage.
%
% EXAMPLE:
% a = 0.95;                  %Coverage set to 955%
% p = 1-a;                   %Coverage rate
% T = 500;                   %Observations
% I = binornd(1,p,T,1);      %Simulates hit-sequence     
% fGMMtest(I,p,1,'uc')       %Calls "GMM-J test" of UC with 1 moment conditions
% fGMMtest(I,p,3,'cc')       %Calls "GMM-J test" of CC with 3 moment conditions
% fGMMtest(I,p,3,'ind')      %Calls "GMM-J test" of Ind with 3 moment conditions
%  
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function [Test, asymptotics, name, varargout] = fGMMtest(I,p,k,criteria,sign,bootstrap)
%Default significance level is set to 5% for the asymptotic critical value
if exist('sign','var') == 0
    sign = 0.05;
end

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

%checks what criteria to test and input
if strcmp(criteria,'ind')
    p = mean(I);
elseif strcmp(criteria,'uc') 
    if k~=1
        disp('Moment restrictions changed to 1 to test unconditional coverage.');
        k=1;
    end
elseif strcmp(criteria,'cc')
    
else
    error('Incorrect input, you must specify either cc, ind or cc criteria');
end

%Calculats asymptotic critical value
% Degrees of freedom for asymptotic distribution is set
if strcmp(criteria,'cc') || strcmp(criteria,'uc')
    df = k;
elseif strcmp(criteria,'ind')
    df = k-1;
end

asymptotics = chi2inv(1-sign,df);
    
%Checks if toolbox is available, then gets durations
if exist('fDurations') == 2
    %Gets Durations
    D = cell2mat(fDurations(I));
else
    error('fDurations not in memory, add the BackTestToolbox before using this test function.');
end

%Checks that significance level is between 0 and 1
if  (sign<=0) || (sign>=1);
    error('Significance level, p, for test is not between 0 and 1.');
end

%checks hit-sequence is of length greater than 2 (a vector)
if (length(I)<2);
    error('Hit-sequence, I, is not of length >1');
end

% %Checks coverage level is a of value between 0 and 1
% if  (p<0) || (p>1);
%     error('Coverage level, p, for test is not between 0 and 1.');
% end

%checks if number of moment conditions is an integer
if (rem(k,1) ~=0)
    error('Number of moment konditions, k, must be a whole number.');
end

if sum(I)>0
    %number of durations
n = length(D);

%Creates the moment conditions matrix
M = OrthPoly(D(:),p,1);
for j=2:k
    M = [M OrthPoly(D(:),p,j)];
end

%sums
Msum = sum(M,1)';

%test statistic is calculated
Test = ((1/sqrt(n))*Msum)'*((1/sqrt(n))*Msum);
elseif sum(I)==0
    disp('Test set to 0, no hits in hit-sequence.');
    Test = 0;
end

%Name of test
name = ['GMM-' criteria '-' num2str(k)];

%Calculates the bootstrapped p-value
if strcmp(bootstrap,'yes')
    varargout{1} = fBootPval(name,Test,I,p,k,criteria,sign,bootstrap);
end

end

%%
% Calculates the Orthonormal polynomials from geometric distrbution
%
% USAGE:
%   M = OrthPoly(d,p,j)
%
% INPUTS:
%   d        -  Duration sequence
%   p        -  Success probability from hit sequence
%   j        -  #Moments used
% 
% OUTPUTS:
%   M       -  Moment condition
%

function M = OrthPoly(D,p,j)

%For j=-1 moment condition is 0
if j == -1; M=0; end

%For j=0 moment condition is 1
if j == 0; M=1; end

%Calculates moment conditions for j larger than 0
if j > 0;
    j=j-1;
    M=(((1-p)*(2*j+1)+p*(j-D+1))/((j+1)*sqrt(1-p))).*OrthPoly(D,p,j)-(j/(j+1))*OrthPoly(D,p,j-1);
end

end


