% Turns the hitsequence input into a sequence indicaitng if each duration is censored, Input may be a
% matrix of hit sequences (with equal length) with 1 series in each column.
% Or it can be a 1 column vector of a single hit-sequence
%
% USAGE:
%   c = fCensoreds(I,D)
%
% INPUTS:
%   I        -  TxN hit-sequence (N sequences of length T)
%   D        -  Duration for each hit-sequence in a cell
%
% OUTPUTS:
%   c        -     Cell of indicators for censoring from the input hit sequence
% 
% EXAMPLE:
% I = binornd(1,p,500,10);  %Gets 10 hit-sequeces of length 500
% d = fDurations(I);        %Gets durations for each sequence
% c = fCensoreds(I,d);      %Gets censoring indicator for each duration
%
% Author:       Thor P. Nielsen (Thor.Nielsen@econ.ku.dk)
% Date:         04-03-2014
% Version:      1.0
%
%%


function c = fCensoreds(I,D)
%memory for censoring series
c = cell(size(D,1),1);

%Create censoring series
for i=1:length(c)
    c{i} = zeros(size(D{i},2),1);
    if I(1,i)==0; c{i}(1) = 1; end       %Checks if first observation is 1
    if I(end,i) ==0; c{i}(end) = 1; end  %Checks if last observation is 1
end

end