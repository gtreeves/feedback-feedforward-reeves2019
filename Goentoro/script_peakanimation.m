% script_peakanimation
%
% This script runs a series of dde simulations of the Goentoro model, and
% plots them, pausing at each step, so that we can visualize how the
% simulations change with a changing parameter.

clear
close all

nK = 200;
K11 = logspace(-2,6,nK)';
K2 = 1e-1;
F0 = 1;
F1 = 10;
n = 1;

y0 = F0./(K11 + F0);
y1 = F1./(K11 + F1);
RHS = y1' - y0';
LHS = (1 + y0'/K2)/F0 - (1 + y1'/K2)/F1;
K3_PA = RHS./LHS;

tauz = 1; tauw = 1; tauy = 1;
K4 = Inf; K8 = Inf;
p = [tauz 0 K2 0 K4 K8 tauw tauy];
theta = [0.5 0.5 0.5];

yesplot = false;
Soln = cell(nK,1);
maxz = 0;
minz = 1;
maxz0 = 0;
minz0 = 1;
for i = 1:nK

	K1 = K11(i);
	K3 = K3_PA(i);
	p(2) = K1;
	p(4) = K3;
	Soln{i} = ftn_rundde(p,theta,F0,F1,n,yesplot);
	z = Soln{i}.y(2,:)';
	maxz = max([z;maxz]);
	minz = min([z;minz]);
	maxz0 = max([z/z(1);maxz0]);
	minz0 = min([z/z(1);minz0]);
end

figure('pos',1.0e+03 *[0.1058    0.2586    1.1264    0.4200])
for i = 1:nK
	
	soln = Soln{i};
	t = soln.x';
	z = soln.y(2,:)';
	
	subplot(1,2,1)
	plot(t,z)
	ylim([0.9*minz 1.1*maxz])
	
	subplot(1,2,2)
	plot(t,z/z(1))
	ylim([0.9*minz0 1.1*maxz0])
	
	pause(0.01)
	F(i) = getframe;
	
	imwrite(F(i).cdata,['Figs_getframe/peakanimation',num2strDU(i,3),'.jpg'])
	
end