% script_Fig1
%
% This script is designed to run all simulations, etc, to get Fig 1 of
% Reeves 2019 (FF/FB).

clear
close all
options = ddeset('RelTol',1e-6);
options2 = optimoptions('lsqcurvefit','display','off');

%
% Fig 1AB are drawn diagrams. So the first one here is Fig 1C. All commands
% are embedded in "ftn_CHE_FFFB".
%
% ftn_CHE_FFFB

%
% Fig 1D: showing a simulation of the FFL of the Goentoro model.
%
t0 = 0; x0 = 1; x1 = 10;
disthand = @(t)ftn_unitstep(t,t0,x0,x1);
n = 1;
theta = 4*[0.5 0.5 0.5];
K1 = 1; K2 = 0.3; K12 = 0.3; K3 = 0.01; K4 = 0.01;
K5 = K1*K3; K6 = K2*K3; K7 = K12*K3;
r = 4; s = 4; tauy = 4;
p = [r K1 K2 K12 K3 K4 s tauy];
tspan = [-10 105];

%
% Initial conditions
%
y0 = x0.^n./(K1.^n + x0.^n);
A0 = x0./K1;
B0 = 1 + A0 + y0./K2 + x0.*y0./K12;
C0 = 1/K3 + x0./K5 + y0./K6 + x0.*y0./K3./K12;
z0 = -(K4.*B0-A0)/2./(B0+C0) + ...
	sqrt((K4.*B0-A0).^2 + 4*(B0+C0).*A0.*K4)/2./(B0+C0);
w0 = z0/K4./(1 + z0/K4);

%
% Perturbed FF/FB sim
%
ftnhand = @ftn_goentoro_dde;
Y0 = [y0; z0; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
t = linspace(tspan(1),tspan(2),500)';
Y = deval(soln,t);
% t = soln.x';
% Y = soln.y';
z = mDU(Y(2,:)'); z = 7*z + 1;
T1 = 4 + 5*(t > 0);

figure('paperpositionm','aut')
plot(t,T1,'Linewidth',2)
hold on
plot(t,z,'linewidth',2)
set(gca,'fontsize',24,'ytick',[])

xlabel('time')
xlim([-5 105])
ylim([0 10])
print(gcf,'Figs/Fig1D_IFFL.jpg','-djpeg','-r150')
print(gcf,'Figs/Fig1D_IFFL.eps','-depsc')





