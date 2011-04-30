%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: gmm1d.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:40 $
%% Version:   $Revision: 1.1 $
%%=============================================================

function [density] = gmm1d(x, w, u, sigma)

n = length(w);
density = zeros(size(x));

for i=1:n
    density = density + w(i)*dnorm(x, u(i), sigma(i));
end
   