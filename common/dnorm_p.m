function [density] = dnorm_p(x,u,sig)

%x= 150
%u = 0
%sig = 87.9005

p = length(x);  % dimensions

% p = length(u)
% k = p*(p+1)/2 = length(sig);

% construct pxp SPD sigma from  sig

a = diag(sig(1:p));
n = eye(p);
i = 0;
for k =1:p-1
   for j = 1: p-k
        i = i+1;
        n(j,j+k) = sig(p+i);
   end  
end
m_sig = n' * a * n;
dist =    (x-u) * m_sig * m_sig * (x-u)' ;
% where x-u is 1xp and m_sig is of pxp, hence dist is a positive scalar
exponent = exp(-dist/2);
density =  det(m_sig) * exponent /  sqrt(2*pi)^p ;

