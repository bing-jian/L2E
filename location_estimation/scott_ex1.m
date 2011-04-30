%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: scott_ex1.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:44 $
%% Version:   $Revision: 1.1 $
%%=============================================================

% consider a sample of size 400 from N(0,1)

n = 400;
x(:,1) = randn(n,1);

% with up to 25 additional samples from a contamination density, N(5,1)

k = [0 20 40 80 120 200];
per = {'0'
       '5%'
       '10%' 
       '20%' 
       '30%' 
       '50%'};
for i=2:6
    x(:,i) = randn(n,1);
    x(1:k(i),i) = 5+randn(k(i),1);
end

mu = mean(x,1);

cellcolors = {
    'k-'
    'r--'
    'cx'
    'm:'
    'g+:'
    'b-.'
};

close all
hold on
tt = linspace(-6,8,120);
for i=1:6
  for t = 1:120
     L2E(t,i) = 1/(2*sqrt(pi)) - 2*mean(dnorm(x(:,i), tt(t), 1));
  end
  plot(tt, 10*L2E(:,i), char(cellcolors(i)), 'linewidth',2);
end
axis([-6,8,-3,3])
hold off
           legend(char(per(1)), char(per(2)), char(per(3)), char(per(4)), char(per(5)), char(per(6)), 'percentage of outliers', 0) 
      xlabel('\mu','fontsize',18, 'fontweight','bold'); 
      ylabel('L2E', 'fontsize',16, 'fontweight','bold');
      %title('L2E', 'fontsize',16, 'fontweight','bold'))
  set(gca,'fontsize',16);

tt = linspace(-10,10,200);
figure
axis([-10,10,-30,0])
hold on
for i=1:6
  plot (tt,  1*(-((mu(i)-tt).^2)/2 - log(2*pi)), char(cellcolors(i)), 'linewidth',2 );
  plot ([mu(i) mu(i)], [-30, -log(2*pi)], 'k-.');
end

hold off
           legend(char(per(1)), char(per(2)), char(per(3)), char(per(4)), char(per(5)), char(per(6)), 'percentage of outliers', 0) 
      xlabel('\mu','fontsize',18, 'fontweight','bold'); 
      ylabel('Log-likelihood', 'fontsize',16, 'fontweight','bold');
      %ylabel('MLE', 'fontsize',16, 'fontweight','bold'))
  set(gca,'fontsize',16);


