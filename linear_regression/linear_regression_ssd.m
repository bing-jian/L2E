% L2E_LINEAR_REGRESS
%
% Author:  Bing Jian (bjian@cise.ufl.edu)
% Date:    May 23, 2004
%
%
% Partial L2E multivariate linear regression
% Purpose:  fit model  Y = a + X*b + ei  by L2E
% the error distribution can be either  N(0,sig^2)  or  w N(0,sig^2)
%
% Inputs:
%           x - n x p input data matrix (p=1 vector OK)
%           y - n responses
%  init_param - input guesses for (a,b,sig) length 1+p+1
% 		defaults = c(mean(y),0,..,0,sd(y))
%       opt_w - False (default e ~ N(0,sig^2) or  True (e ~ w N(0,sig^2))
%      init_w - fixed value for w  or initial guess for w (see opt_w)
%
%
% Output:      list containing estimated coefficients
%          $a - intercept
%          $b - vector of estimated beta's
%        $sig - estimated standard deviation of residuals
%          $w - estimated weight (note:  w can be greater than 1)
%        $res - residual vector
%
% See also l2e_linear_regress_criterion
%
%
% Reference:
%   [1] Scott, D.W. (2001), Parametric Statistical Modeling by Minimum
%    Intergrated Square Error, Technometrics, 43, 274-285
%   [2] S-Plus code by Scott, D.W.,
%     http://www.stat.rice.edu/people/scottdw/public_html/code/l2e/
%


function [a, b, res] = ssd_linear_regress( X, Y , init_param)

%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: linear_regression_ssd.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:44 $
%% Version:   $Revision: 1.1 $
%%=============================================================


%% check input
n = size(X,1);  % size of sample
p = size(X,2);  % number of variables

if (( n==1) && (p>1)) % consider X as a row vector
    X = X';
    n = p;
    p = 1 ;
    Y = Y';
end

if length(Y) ~= n
    disp('X and Y differ on n');
    return
end



if nargin < 3
	% set default initial guess of parameters
	a = mean(Y);
	b = zeros(p,1);
	init_param = [a;b];   % construct column vector of initial parameters
end

if length(init_param) ~= (p+1)
        disp('length of init_param is incorrect');
	return
end


options=optimset('display','final', 'MaxFunEvals', 1000000);
[final_param, min_value] = fminunc(@ssd_linear_regress_criterion, init_param, options, X, Y);

a = final_param(1);
b = final_param(2:p+1);
res = Y - a - X*b;

