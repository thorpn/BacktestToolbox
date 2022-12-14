% Calculates bootstrapped p-value using the method of Dufour (2006)
% 
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     04-06-2014
% Version:  1.0
%
%%

function Pval = fBootPval(varargin)

name = varargin{1};         %name of test
test = varargin{2};         %test value
I = varargin{3};            %Hit-seq
p = varargin{4};            %coverage rate
n = 2000;                   %standard number of bootstraps is set to 2000
T = size(I,1);              %Number of observations

%Gets test name
TestFunction = fFindTestFunction(name);
TestFunction = str2func(TestFunction);

%Generates samples under the null
I = binornd(1,p,T,n); 

TestValues = NaN(n,1);
for i=1:n
    varargin{3} = I(:,i);
    TestValues(i) =  TestFunction(varargin{3:end-1});
end

s0 = test;
u0 = unifrnd(0,1);
u = unifrnd(0,1,n,1);

LessThans = TestValues<s0;
Equalz = abs(TestValues-s0)<0.000000000001;
Ghat = 1 - (1/n)*sum(LessThans+Equalz)+(1/n)*sum(Equalz.*(u0>=u));

Pval = (n*Ghat+1)/(n+1);

end

%%
% Recognizes the test function from a list
%
% USAGE:
%   estFunction = fFindTestFunction(name)
%
% INPUTS:
%   name        -  Name output from one of the backtest functions
%
% OUTPUTS:
%   TestFunction - Function name of backtest function
%
% Author:       Thor Nielsen
% Date:         19/08/2014

function TestFunction = fFindTestFunction(name)

if strcmp(name,'PF');
    TestFunction = 'fPFtest'; 
elseif strcmp(name,'Z-test');       
    TestFunction = 'fztest'; 
elseif  strcmp(name,'TUFF');
    TestFunction = 'fTUFFtest';
elseif strcmp(name,'Joint'); 
    TestFunction = 'fJointtest';
elseif strcmp(name,'DurCtest');
    TestFunction = 'fDurCtest';
elseif strcmp(name(1:7),'Markov-');  
    TestFunction = 'fGeneralizedMarkovtest';
elseif strcmp(name(1:4),'DurM');  
    TestFunction = 'fDurM'; 
elseif  strcmp(name(1:3),'LB-');      
    TestFunction = 'fLBtest'; 
elseif   strcmp(name(1:8),'DurDtest'); 
    TestFunction = 'fDurDtest'; 
elseif  strcmp(name(1:4),'GMM-');     
    TestFunction = 'fGMMtest';
elseif strcmp(name(1:3),'DQ-');      
    TestFunction = 'fDynamicQuantileTest'; 
end

end



