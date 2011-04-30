% L2E_LINEAR_REGRESS_TESY
%
% Script for testing multivariate linear regress using L2E
%
% Author:   Bing Jian
% Date:     May 23, 2004
%
% See also l2e_linear_regress, l2l2e_linear_regress_criterion
%
% Reference:
%   [1] Scott, D.W. (2001), Parametric Statistical Modeling by Minimum
%    Intergrated Square Error, Technometrics, 43, 274-285
%   [2] S-Plus code by Scott, D.W.,
%     http://www.stat.rice.edu/people/scottdw/public_html/code/l2e/

%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: linear_regression_gamma_test.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:43 $
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
%    beta = sig*sig/(a+eps);
    beta = 1;
	init_param = [a;b;alpha];   % construct column vector of initial parameters
        opt_w = false;
        init_w = 1;

[a, b, alpha, w, res] = gamma_linear_regress( X, Y , init_param, opt_w, init_w)
a
b
alpha
beta
w

close all
plot(X,Y, 'o');
hold on
plot(X, a + X*b,'g');
hold on
plot(X, 2 - 3*X, 'r');
hold off
