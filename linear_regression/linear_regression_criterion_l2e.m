% L2E_LINEAR_REGRESS_CRITERION
%
% Author:  Bing Jian (bjian@cise.ufl.edu)
% Date:    May 23, 2004
%
%
% Partial L2E multivariate linear regression
% Purpose:  fit model  Y = a + X*b + ei  by L2E
% the error distribution can be either  N(0,sig^2)  or  w N(0,sig^2)
%
% Given following input arguments:
%       X       --    Input     n*p matrix  where p is number of variables
%       Y       --    Response  n*1 vector  where n is number of samples
%       a, b    --    Coefficients of linear fitting function Y = a + X*b + ei
%                     where "a" is scalar and b is p*1 vector
%       w, sig  --    Parameters of normal distribution for error
%
% The criterion function to be minimized can be expressed as
%       f(w,sig, ei) = w^2/(2*sqrt(pi)*sig) - 2*w*mean(dnorm(ei,0,sig));
%       where ei = Y - a - X*b   and
%             dnorm(x,mu,sig) is pdf of N(mu, sig)
%
% See also l2e_linear_regress
%
% Reference:
%   [1] Scott, D.W. (2001), Parametric Statistical Modeling by Minimum
%    Intergrated Square Error, Technometrics, 43, 274-285
%   [2] S-Plus code by Scott, D.W.,
%     http://www.stat.rice.edu/people/scottdw/public_html/code/l2e/
%

function [val] = l2e_linear_regress_criterion(param, X, Y, opt_w, init_w)
%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: linear_regression_criterion_l2e.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:42 $
%% Version:   $Revision: 1.1 $
%%=============================================================

n = size(X,1);
p = size(X,2);

a = param(1);
b = param(2:p+1);
sig = param(p+2);

if (opt_w)           % if opt_w is on
  w = param(p+3);    % then w is also a parameter to be found
else
  w = init_w;        % otherwise w is fixed to init_w
end


ei = Y - a - X*b;
val = w^2/(2*sqrt(pi)*sig) - 2*w*mean(dnorm(abs(ei),0,sig));
val = 100*val;

