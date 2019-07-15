function[sum] = geometric(r,n)
% Calculates the sum  of the geometric series r^n
% Note that when -1<r<1 then the series converges to 1/(1-r)

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

count = 1; sum = 1;

while count <= n
    sum = sum + r^count;
    count = count + 1;
end