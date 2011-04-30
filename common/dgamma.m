function [density] = dgamma(x,alpha,beta,u)

%x: random variable
%alpha:  shape parameter
%beta: scale parameter: default = 1
%u: location parameter: default = 0;

if nargin<4
    u = 0;
end
if nargin<3
    beta = 1;
end

if (x<u)
    disp('input x should not be less than location parameter u');
    return
end

dist = (x-u)/beta;
density = power(dist,alpha-1).*exp(-dist)./(beta*gamma(alpha));
