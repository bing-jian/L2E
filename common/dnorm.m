function [density] = dnorm(x,u,sig)

%x= 150
%u = 0
%sig = 87.9005
exponent = (-(x-u).^2)/(2*sig*sig);
density = exp(exponent);
density = density / (sig * sqrt(2*pi));