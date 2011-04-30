%%=============================================================
%% Project:   L2E
%% Module:    $RCSfile: plot_gmm1d.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:40 $
%% Version:   $Revision: 1.1 $
%%=============================================================

u = [0,5];
sigma = [1,1];
x = -5:0.1:10;

per = {'0'
       '5%'
       '10%' 
       '20%' 
       '30%' 
       '50%'};

w1 = 1 - [0, 0.05, 0.1, 0.2, 0.3, 0.5];
w2 = 1 - w1;

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
   w = [w1(i), w2(i)];
   density = gmm1d(x,w, u,sigma);
   plot(x,density,char(cellcolors(i)), 'linewidth',2);
end
hold off
      legend(char(per(1)), char(per(2)), char(per(3)), char(per(4)), char(per(5)), char(per(6)), 'percentage of outliers', 0) 
      xlabel('x','fontsize',18, 'fontweight','bold'); 
      ylabel('f(x)', 'fontsize',16, 'fontweight','bold');
      %title('L2E', 'fontsize',16, 'fontweight','bold'))
  set(gca,'fontsize',16);



