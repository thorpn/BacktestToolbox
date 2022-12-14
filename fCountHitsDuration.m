% Counts N00, N01, N10, N11. This code was used to generate the MEX file
% which is much faster. It is called from fDurM.m
%
% USAGE:
%   [N11, N10, N00, N01]  = fCountHitsDuration(I,states)
%
% INPUTS:
%   I         -  Hit-sequence, I, column vector
%   states    -  Number of lags against which we want power
%
% OUTPUTS:
%   [N11, N10, N00, N01] -  Counting variables of hit combinations
%
% Author:   Thor P. Nielsen (econ.ku.dk/pajhede)
% E-mail:   thorpn86@gmail.com
% Date:     14-08-2014
% Version:  1.0
%%

function [N11, N10, N00, N01]  = fCountHitsDuration(I,states)
[N11, N10] = deal(zeros(states,1));
N00 = 0;    N01 = 0; 

for i=(states+1):(length(I))

    %%Picks state
    Limes = NaN(states,1);
    Limes(:) = 1:states;

    for j=1:states
        Interval     = i-Limes(j);
        PrevInterval = Interval+1;

        %First states doesnt condition
        if (j==1 && I(Interval)==1)
                N10(j)=N10(j)+(1-I(i));
                N11(j)=N11(j)+I(i);

        %Other states condition on newers observations being 0
        elseif (sum(I(PrevInterval:(i-1)))==0 && j~=1 && I(Interval)==1);
                N10(j)=N10(j)+(1-I(i));
                N11(j)=N11(j)+I(i);
        end
    end

    %estimates based on no hits previously
    if sum(I((i-states):(i-1))) == 0
         N00=N00+(1-I(i));
         N01=N01+I(i);
    end
end

end


