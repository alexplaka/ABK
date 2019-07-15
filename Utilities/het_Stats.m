function s_het = het_Stats(PHM,flag)
% indHet calculates the index of heterogeneity, psi, given the
% Population Heterogeneity Matrix (PHM). The PHM can be for a
% 1st, or 2nd order process. The argument "flag" should be 
% set to 1 (or a nonzero value) for homologous 2nd order processes.
% In that case, the PHM is a square matrix whose diagonal values
% are set to 0 (if they haven't already been set to 0). This is 
% because an agent can not intercat with itself in such a process.
% The function returns the structure s_het with fields:
% 1) 'k_mean': the mean of nonzero (or non-NaN) k values in the PHM. 
% 2) 'k_sdev': the standard deviation of nonzero (or non-NaN) k values in the PHM. 
% 3) 'psi': the index of heterogeneity, psi.
% 4) 'chi': the subspecies fractional abundances vector, chi.
% 5) 'S': information-theoretic measure of population heterogeneity
% 6) 'k_unique': a vector with the corresponding unique k
%       values in the PHM.

% TODO: May need to make adjustments to the following code 
% wrt assignment of 0's or NaN's in PHM (esp. for homologous 2nd order processes).

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

if nargin == 1
    flag = 0;
end

if flag ~= 0
    for w=1:size(PHM,1)
       PHM(w,w) = 0; 
    end   
end

[k, ip, ic] = unique(PHM);      % Find unique k values in PHM

% Remove k=0 values (using 0 instead of NaN because "unique" function 
% considers NaN entries as distinct.)
if k(1) == 0                    
    k = k(2:end);
    ip = ip(2:end);
    ic = ic(2:end);
end

C = size(k,2);                  % Species richness (ie, how many subspecies are there?)

n = accumarray(ic,1);           % n is a vector containing the subspecies population sizes

N = sum(n);                     % Total number of agents

chi = n / N;                    % Vector of subspecies fractional abundances


l = 0;
s = 0;

for i=1:C
    l = l + n(i) * (n(i)-1);
    s = s - chi(i) * log(chi(i));  % Information-theoretic measure of pop. heterogeneity
end

L = l / (N * (N-1));

psi = 1 - L;

k_unique = k';


if sum(size(PHM)==1) ~= 0   % if the PHM is for a 1st order process (ie, is an array or vector)
    k_mean = mean(PHM);
    k_sdev = std(PHM);
else                        % else PHM is for a 2nd order process
    PHM(PHM==0) = NaN;      % replace 0's with NaN's
    k_mean = mean(mean(PHM,'omitNaN'),'omitNaN');
    k_sdev = std(std(PHM,'omitNaN'),'omitNaN');
end

s_het = struct('k_mean',k_mean,'k_sdev',k_sdev,'psi',psi,...
    'S',s,'chi',chi,'k_unique',k_unique);