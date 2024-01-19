function y = negFb2_dde_hist(t)

% Author: Alex Plakantonakis,   Copyright (c) 2019-2024.           License: GNU GPLv3

global DEhistory;

y = deval(DEhistory,t);