function[a] = alex_factorial(n)
% This function calculates the factorial of any number n.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

if n==0
    a=1;
end

count = 1;
a = 1;

while count < n
    a = a * (count + 1);
    count = count + 1;
end
