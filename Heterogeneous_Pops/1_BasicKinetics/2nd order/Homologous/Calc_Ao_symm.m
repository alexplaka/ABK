% The following applies to the process 2A --> X for heterogeneous populations
% of A, where the PHM has size Ao x Ao, where Ao is the initial 
% population of A. Because k_ii for i=1:Ao (the diagonal entries) don't
% make sense to consider in the ABK implementation (an agent can not
% interact with itself in a second order process for a population
% with a discrete number of agents) the number of available 
% interactions is thus Ao*(Ao-1).

% This script calculates the (initial) population size, such that,
% given a specified number of subinteraction groups s, there is an 
% equal number of (reciprocal) agent-agent pairings within each
% distinct subinteraction group. In other words, the reported choices 
% of Ao guarantee a symmetrically-distributed heterogeneous population.

% The relationship that provides the caculation of Ao is
%               Ao*(Ao-1) / (2*s) = 2*k
% where k=0,1,2,... (whole numbers). Evaluation of the left side should 
% yield an even number for this intrapopulation symmetry to exist. 
% The resulting quadratic equation is what the code below solves for
% using integer values of k.The division by 2 on the left-hand side 
% assumes that agent-agent interactions are reciprocal, as represented 
% in the PHM (e.g., k_12 = k_21).

% Author: Alex Plakantonakis,   Copyright (c) 2019.           License: GNU GPLv3

clear;   clc;

s = 3;
disp(['Number of subinteraction groups = ' num2str(s)]);

for k=0:100
    Ao = (1 + sqrt(1 + 16*s*k)) / 2;
    if Ao == floor(Ao)            % equivalently, exp(2*pi*i*N) == 1 only if N is an integer 
        disp(['k=' num2str(k) ' ,    Ao = ' num2str(Ao)]);
    end
end        


% Note: 
%  - For s = 2, the computed Ao values are 
% Ao = 1, 8,9, 16,17, 25,25, 32,33, ...