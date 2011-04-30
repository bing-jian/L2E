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

function [val] = l2e_robust_mean_criterion(params, X, opt_w, init_w)
%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: l2e_robust_mean_criterion.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:44 $
%% Version:   $Revision: 1.1 $
%%=============================================================

if (nargin<4)
  init_w = 1;
end
if (nargin<3)
 opt_w = false;
end
%params
n = size(X,1);  % number of points
p = size(X,2);  % dimensions

k = p*(p+1)/2;
center = params(1:p); % 1xp
%mu = param(p+1:2*p);
% let's fix mu to 0
mu = zeros(1,p);
sig =  params(p+1:p+k);
% construct pxp SPD sigma from sig
% such that   sig*sig*SIGMA = I
a = diag(sig(1:p));
% nn = eye(p);
% i = 0;
% for k =1:p-1
%    for j = 1: p-k
%         i = i+1;
%         nn(j,j+k) = sig(p+i);
%    end
% end
% m_sig = nn' * a * nn;

m_sig = a;
if (opt_w)           % if opt_w is on
  w = params(p+k+1);    % then w is also a parameter to be found
else
  w = init_w;        % otherwise w is fixed to init_w
end

ei = X - ones(n,1)*center;  % nxp

val = w^2/(2*sqrt(pi))^p;

val = val * det(m_sig);

sum = 0;
for i=1:n
  sum = sum + dnorm_p(ei(i,1:p), mu, sig);
end

val = val - 2*w*sum/n;
val = 10000*val;

