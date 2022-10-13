% script_FigS1
%
% This script is designed to run the simulation from Fig S1 of Reeves 2019
% (FF/FB).
%
% This is simply a simulation of the system in Region 1.

clear
close all
options = ddeset('RelTol',1e-6);

t0 = 0; F0 = 1; F1 = 10;
n = 1;
theta = [0.5 0.5 0.5];
K1 = 0.01; K2 = 100; K12 = 0.3; K3 = Inf; K4 = Inf;
r = 1; s = 1; tauy = 1;
p = [r K1 K2 K12 K3 K4 s tauy];
tspan = [-5 10];
yesplot = false;
soln = ftn_rundde(p,theta,F0,F1,n,yesplot,tspan);
disthand = @(t)ftn_unitstep(t,t0,F0,F1);
t = soln.x';
Y = soln.y;
y = Y(1,:)';
z = Y(2,:)';
x = disthand(t);

figure('paperpositionm','aut')
% plot(t,x,'Linewidth',2)
plot(t,y,'linewidth',2)
hold on
plot(t,z,'linewidth',2)

set(gca,'fontsize',24)

legend('y','z')
% set(gca,'ytick',0.94:0.01:1)
print(gcf,'Figs/FigS1.eps','-depsc')
print(gcf,'Figs/FigS1.jpg','-djpeg','-r150')

ylim([0 1])
set(gca,'ytick',0:0.2:1)
print(gcf,'Figs/FigS1_zoomout.eps','-depsc')




