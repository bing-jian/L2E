function r=gammaplot(xmin,xmax,a,b,col)
%gammaplot(xmin,xmax,a,b,col)
%plots a GAMMA distribution with parameter a and b > 0, in the interval [xmin,xmax]
%col is a text string which indicates graph colour (default col='b')
%Copyright gianlucabaio2002

if nargin==4
   col='b';
end

x=[xmin:.1:xmax];
const1=b^a;
const2=gamma(a);

sum = 0;
for i=1:length(x)
   y(i)=const1/const2*x(i)^(a-1)*exp(-x(i)*b);
   z(i) = dgamma(x(i),a,b);
   sum = sum+z(i);
end


plot(x,y,col)
hold on
plot(x,0.2+z,'g')
hold off

r = sum;