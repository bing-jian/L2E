% L2E_LINEAR_REGRESS_TEST
%
% Script for testing multivariate linear regress using L2E
%
% Author:   Bing Jian
% Date:     May 23, 2004
%
% See also l2e_linear_regress, l2l2e_linear_regress_criterion
%
% Reference: 
%   [1] Scott, D.W. (2001), Parametric Statistical Modeling by Minimum
%    Intergrated Square Error, Technometrics, 43, 274-285
%   [2] S-Plus code by Scott, D.W., 
%     http://www.stat.rice.edu/people/scottdw/public_html/code/l2e/ 

%%=====================================================================
%% Project:   L2E
%% Module:    $RCSfile: linear_regression_l2e_test.m,v $
%% Language:  MATLAB
%% Author:    $Author: bjian $
%% Date:      $Date: 2008/12/09 22:52:43 $
%% Version:   $Revision: 1.1 $
%%=====================================================================


clear
X = linspace(-10,10,100);
Y = 2 - 3*X;

n = length(Y);
noise = 2*(rand(1,n) - 0.5);        % noise
Y = Y + noise;

% outlier
Y(26:30) = 30 + 2*rand(1,5);  
Y(61:70) = 20 + 3*rand(1,10);
Y(81:90) = 40 + 3*rand(1,10);
Y(6:10) = 0 + 3*rand(1,5);
Y(96:100) = 0 + 3*rand(1,5);


p = 1;
a = mean(Y);
b = zeros(p,1);
sig = std(Y);
init_param = [a;b;sig];   %
opt_w = true;
init_w = 1;
          
[a, b, sig, w, res1] = l2e_linear_regress( X, Y , init_param, opt_w, init_w)
a
b
sig
w

close all
export_fig_path = './';

set(gcf, 'DefaultTextFontWeight','bold','DefaultTextFontSize',24, 'DefaultTextInterpreter', 'latex') 
set(gcf, 'DefaultAxesFontSize', 24);

plot(X,Y, 'k*');

hold on
%plot(X, 2 - 3*X, 'b:', 'linewidth', 3);
plot(X, a + X*b,'r', 'linewidth', 3);

[a, b, res2] = ssd_linear_regress( X, Y )
plot(X, a + X*b,'g','linewidth', 3);
legend('Data','L_2E','SSD');
hold off
axis off
filename = 'l2e_vs_ssd';
eval(['print -dpng ',sprintf('%s',export_fig_path),filename,'.png']);
eval(['print -deps2 ',sprintf('%s',export_fig_path),filename,'_gray.eps']);
eval(['print -depsc2 ',sprintf('%s',export_fig_path),filename,'_color.eps']);

figure

set(gcf, 'DefaultTextFontWeight','bold','DefaultTextFontSize',20, 'DefaultTextInterpreter', 'latex') 
set(gcf, 'DefaultAxesFontSize', 24);

bin1 = hist(res1,100);
bar(bin1/100);
hold on
t = linspace(min(res1),max(res1),100);
plot(w*dnorm(t,0,sig),'r','linewidth',2)
legend('histogram of residues','estimated partial Gaussian');
axis([0,120,0,0.4])
hold off
filename = 'fitting_residue';
eval(['print -dpng ',sprintf('%s',export_fig_path),filename,'.png']);
eval(['print -deps2 ',sprintf('%s',export_fig_path),filename,'_gray.eps']);
eval(['print -depsc2 ',sprintf('%s',export_fig_path),filename,'_color.eps']);

