%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: nnn.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:44 $
%% Version:   $Revision: 1.1 $
%%=============================================================

clear
X = linspace(-10,10,100);
Y = 2 - 3*X;

n = length(Y);
noise = 2*(rand(1,n) - 0.5);        % noise
Y = Y + noise;

% outlier
Y(21:30) = 30 + 2*rand(1,10);
Y(61:70) = 20 + 3*rand(1,10);
Y(81:90) = 40 + 3*rand(1,10);
Y(6:10) = 0 + 3*rand(1,5);
Y(96:100) = 0 + 3*rand(1,5);


        p = 1;
	a = mean(Y);
	b = zeros(p,1);
	sig = std(Y);
    alpha = (a/sig)^2;
    beta = sig*sig/(a+eps);
    beta = 1;
	init_param = [a;b;alpha;beta];   % construct column vector of initial parameters
        opt_w = false;
        init_w = 1;
