% script_FigS3
%
% In this script, we will run several simulations of FF/FB with n = 2, K1 =
% 1, K2 = 0.1, and show that the peak value is independent of K3,K4.


clear
close all
options = ddeset('RelTol',1e-6);

t0 = 0; F0 = 1; F1 = 10;
n = 2;
theta = [0.5 0.5 0.5];
K1 = 1; K2 = 0.1; 

y0 = F0.^n./(K1.^n + F0.^n);
y1 = F1.^n./(K1.^n + F1.^n);
RHS = y1'.^n - y0'.^n;
LHS = (1 + ((1./K2)*y0').^n)/F0.^n - (1 + ((1./K2)*y1').^n)/F1.^n;
% K3star = (repmat(RHS,nK2,1)./LHS).^(1/n);
K12 = (RHS./LHS).^(1/n);

K3 = Inf; K4 = Inf;
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 15];
yesplot = false;


K33 = [0.01 0.1 1];
K44 = [0.01 0.1 1];
T = cell(length(K33),length(K44));
Z = cell(length(K33),length(K44));
W = cell(length(K33),length(K44));
figure('paperpositionm','aut')
for i = 1:length(K33)
	for j = 1:length(K44)
		K3 = K33(i);
		K4 = K44(j);
		
		p = [tauz K1 K2 K12 K3 K4 tauw tauy];
		soln = ftn_rundde(p,theta,F0,F1,n,yesplot,tspan);
		disthand = @(t)ftn_unitstep(t,t0,F0,F1);
		t = soln.x';
		T{i,j} = t;
		Y = soln.y;
		z = Y(2,:)';
		Z{i,j} = z;
		w = Y(3,:)';
		W{i,j} = w;
		
		plot(t,z/z(1))
		hold on
	end
end
set(gca,'fontsize',24)

% set(gca,'ytick',0.94:0.01:1)
print(gcf,'Figs/FigS3_Z.eps','-depsc')
print(gcf,'Figs/FigS3_Z.jpg','-djpeg','-r150')




figure('paperpositionm','aut')
for i = 1:length(K33)
	for j = 1:length(K44)
		
		plot(T{i,j},W{i,j}/W{i,j}(1))
		hold on
	end
end

set(gca,'fontsize',24)

% set(gca,'ytick',0.94:0.01:1)
print(gcf,'Figs/FigS3_W.eps','-depsc')
print(gcf,'Figs/FigS3_W.jpg','-djpeg','-r150')


