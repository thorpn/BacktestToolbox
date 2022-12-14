% Counts n00, n01, n10, n11. This code was used to generate the MEX file
% which is much faster. It is called from fGeneralizedMarkovtest.m
%
% USAGE:
%   [n00, n01, n10, n11]  = fCountHits(I,states)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   states    -  Number of lags against which we want power
%
% OUTPUTS:
%   [n00, n01, n10, n11] -  Counting variables of hit combinations
%
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     14-08-2014
% Version:  1.0
%%

function [n00, n01, n10, n11]  = fCountHitsGeneralized(I,states)

[n00, n01, n10, n11] = deal(0);

for i=(states+1):(length(I))

    if sum(I((i-states):(i-1))) == 0
         n00=n00+(1-I(i));
         n01=n01+I(i);
    end

    if sum(I((i-states):(i-1))) > 0
         n10=n10+(1-I(i));
         n11=n11+I(i);
    end

end

end


