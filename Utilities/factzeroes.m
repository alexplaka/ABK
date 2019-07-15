function [z] = factzeroes(n)
% This function calculates the number of terminal zeros in the result of n!
% The value of the factorial is not calculated. Instead, the total
% instances of the prime factors 2 and 5 are counted and used to calculate
% their product, whose number of terminal zeros is returned. 

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

array = 1:n;
t = 0;     f = 0;

for x=1:size(array,2)
    twos = find(factor(array(x))==2);
    t = t + size(twos,2);
    fives = find(factor(array(x))==5);
    f = f + size(fives,2);
end

factors25 = 2^t * 5^f;
str = num2str(factors25);
terminalz = regexp(str,'0+$');
if isempty(terminalz) == 1      
    z = 0;
else
    z = size(str,2) - terminalz + 1;
end
