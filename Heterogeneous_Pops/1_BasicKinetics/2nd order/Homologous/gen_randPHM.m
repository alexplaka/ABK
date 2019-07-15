function [k , abundance] = gen_randPHM(Ao,k_sub)
% Given the rate constants for the subinteraction groups,
% this function generates a PHM with the rate constants for 
% randomly chosen interacting agent-agent pairs. The interactions
% are assumed to be reciprocal: k(i,j) = k(j,i)

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

s = size(k_sub,2);              % number of subinteraction groups
k = zeros(Ao);

for i=1:Ao
    for j=i+1:Ao
           k(i,j) = k_sub(ceil(s*rand));
           k(j,i) = k(i,j);
    end
end

for x=1:s
    abundance(x) = numel(find(k==k_sub(x))) / (Ao*(Ao-1));
    disp(['SubInt-' num2str(x) ' Abundance: ' num2str(abundance(x))]);
end
