% Turns the hit-sequence input (I) into a sequence of durations. Input may be a
% matrix of hit sequences (with equal length) with 1 series in each column.
% Or it can be a 1 column vector of a single hit-sequence
%
% USAGE:
%   d = fDurations(I)
%
% INPUTS:
%   I        -  TxN hit-sequence (N sequences of length T)
%
% OUTPUTS:
%   d     -     n Cells of durations for each series
% 
% EXAMPLE:
% I = binornd(1,p,500,10);  %Gets 10 hit-sequeces of length 500
% d = fDurations(I);        %Gets durations for each sequence
% 
% Author:       Thor P. Nielsen (Thor.Nielsen@econ.ku.dk)
% Date:         04-03-2014
% Version:      1.0
%
%%

function d = fDurations(I)
[Observations, Series] = size(I);

%Memory for durations
d = cell(size(I,2),1);

%loops over each series
for k=1:Series
    indices =zeros(1,1);
    
    %loops over observations in series, if hit is found it is added to
    %indices of hits
    for l=1:Observations
       if I(l,k)==1; indices = [indices l]; end
    end
    
    %If first observation is hit, discard 0 as first observation
    if I(1,k)==1; indices=indices(2:end); end;
    
    %If last observation is hit, add last observation to indices
    if I(Observations,k)==0; indices=[indices Observations]; end;
    
    %Take differences to get durations
    d{k} = diff(indices);
    
end
    
end
