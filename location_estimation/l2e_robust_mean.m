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


function [final_param, min_value] = l2e_robust_mean(X, init_param, opt_w, init_w)
%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: l2e_robust_mean.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:44 $
%% Version:   $Revision: 1.1 $
%%=============================================================

%% check input
% n = size(X,1);  % size of sample
p = size(X,2);  % number of variables

if nargin<4
       init_w = 1;
end

if nargin < 3
       opt_w = false;
end

if nargin < 2
	% set default initial guess of parameters
	center = mean(X);
	%center = 100*ones(1,p);
	mu = zeros(1,p);
	%m_sig = sqrt(inv(cov(X)));
	k = p*(p-1)/2;
	%sig = [svd(m_sig)' zeros(1,k)]
	sig = [ones(1,p) zeros(1,k)];
	%center = [100 100]
	init_param = [center sig] ;  % construct column vector of initial parameters
end

k = p*(p+1)/2;
if (length(init_param) ~= (p+k))
        disp('length of init_param is incorrect');
	return
end

if opt_w
     init_param = [init_param init_w];
end



%options=optimset('display','final', 'MaxFunEvals', 1000000);
options=optimset('display','Iter', 'MaxFunEvals', 1000, 'TolFun',1e-10, 'TolX',1e-10);

pdim = length(init_param);
Lb = -Inf*ones(1,pdim);
Ub = Inf*ones(1,pdim);
Lb(p+1:p+p) = exp(-5);  % sigma >0 %% this constraint on sigma is very important !!
Ub(p+1:p+p) = exp(5);

[final_param, min_value] = fmincon(@l2e_robust_mean_criterion, init_param, [],[], [],[], Lb, Ub, [], options, X,  opt_w, init_w);

init_param(1:p)

final_param
final_param(1:p)
final_param
%mu = final_param(p+1:2*p)

sig = final_param(p+1:p+k)
if opt_w
  w = final_param(p+k+1)
else
  w = init_w
end

min_value

