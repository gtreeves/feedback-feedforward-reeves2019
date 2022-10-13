% script_Fig4
%
% This script is designed to run all simulations, etc, to get Fig 4 of
% Reeves 2019 (FF/FB).
%
% What do we want in Fig 4? This fig is for a FB system as baseline.
%
% Part (a): show how FB alone results in offset, and compare it to the
% FF/FB system.  Need to figure out how to choose K4,K8 and K2,K3 for this
% figure.
%
% Part (b): Show an unstable system with the oscillatory dynamics for the
% FB only system, and also how adding FF stabilizes.
%
% Part (c): Maybe show the same plane and insets for a combined FF/FB
% system. In which case, how do we choose K2? K3? Still let K1 = 1.
%
% If the above gets reserved for the supplement, then what's below will
% actually be Part (c) instead of (d).
%
% Part (d): Now show the K4K8 plane with Alpha1 - alpha2, which also shows
% that for some values of K4,K8, you have unstable systems. 
%
% Maybe hold the Adaptation, Peak, Other peak, and SNR pcolor contours for
% Fig 5?

clear
close all
options = ddeset('RelTol',1e-6);


%% ========================================================================
% First, we will simulate a FB only system with delay and n = 2.
% =========================================================================
%{

t0 = 0;
F0 = 1;
F1 = 10;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Base param values
%
n = 2;
tauz = 1; tauw = 1; tauy = 1;
theta = [0.5 0.5 0.5];
yesplot = true;
tspan = [0 30];
K1 = 1; K2 = 0.1; K3 = 0.3; K4 = 0.3;

y0 = F0.^n./(K1.^n + F0.^n);
y1 = F1.^n./(K1.^n + F1.^n);
RHS = y1'.^n - y0'.^n;
LHS = (1 + ((1./K2)*y0').^n)/F0.^n - (1 + ((1./K2)*y1').^n)/F1.^n;
% K3star = (repmat(RHS,nK2,1)./LHS).^(1/n);
K12 = (RHS./LHS).^(1/n);

p = [tauz K1 K2 K12 K3 K4 tauw tauy];
p1 = p;
p2 = p1; p2(5) = Inf; p2(6) = Inf;

%
% FB only param values
%
figure
p(3) = Inf; p(4) = Inf;
soln = ftn_rundde(p,theta,F0,F1,n,yesplot,tspan);
set(gcf,'paperpositionmode','auto')
set(gca,'fontsize',24)

%
% FF/FB param values
%
yesplot = false;
soln1 = ftn_rundde(p1,theta,F0,F1,n,yesplot,tspan);
hold on
plot(soln1.x,soln1.y(2,:))%/soln1.y(2,1)*soln.y(2,1))

% %
% % FF/FB param values
% %
% yesplot = false;
% soln2 = ftn_rundde(p2,theta,F0,F1,n,yesplot,tspan);
% hold on
% plot(soln2.x,soln2.y(2,:))%/soln1.y(2,1)*soln.y(2,1))
print(gcf,'Figs/negative feedback simulation.eps','-depsc')
print(gcf,'Figs/negative feedback simulation.jpg','-djpeg','-r150')

%}









%% ========================================================================
% Now we will simulate FB only in an unstable regions, & compare to FFFB
% =========================================================================
%{

t0 = 0;
F0 = 1;
F1 = 10;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Base param values
%
n = 2;
tauz = 1; tauw = 1; tauy = 1;
theta = [0.5 0.5 0.5];
yesplot = true;
tspan = [0 30];
K1 = 1; K2 = 0.1; K3 = 0.03; K4 = 0.3;

y0 = F0.^n./(K1.^n + F0.^n);
y1 = F1.^n./(K1.^n + F1.^n);
RHS = y1'.^n - y0'.^n;
LHS = (1 + ((1./K2)*y0').^n)/F0.^n - (1 + ((1./K2)*y1').^n)/F1.^n;
% K3star = (repmat(RHS,nK2,1)./LHS).^(1/n);
K12 = (RHS./LHS).^(1/n);

p = [tauz K1 K2 K12 K3 K4 tauw tauy];
p1 = p;

%
% Unstable simulation
%
figure
p(3) = Inf; p(4) = Inf;
soln = ftn_rundde(p,theta,F0,F1,n,yesplot,tspan);
set(gcf,'paperpositionmode','auto')
set(gca,'fontsize',24)



%
% FF/FB param values
%
yesplot = false;
soln1 = ftn_rundde(p1,theta,F0,F1,n,yesplot,tspan);
hold on
plot(soln1.x,soln1.y(2,:))%/soln1.y(2,1)*soln.y(2,1))
ylim([0 0.21])
print(gcf,'Figs/negative feedback simulation unstable.eps','-depsc')
print(gcf,'Figs/negative feedback simulation unstable.jpg','-djpeg','-r150')

%}


%% ========================================================================
% Now we will do a pcolor_contour of the principal eigenvalue of FB only,
% vs K4,K8, then FF/FB, then the difference.
% =========================================================================
% {

clear
load Mat/2019-01-04_ss_K3K4sweep
load Mat/2018-11-30_eigval_of_J

%
% Now that we have J_out, we can convert those numbers into eigenvalues.
%
Alpha1 = interp1(J,Alpha,J_out(:,:,1));
Alpha2 = interp1(J,Alpha,J_out(:,:,2));

%
% FB only
%
figure
pcolor_contour(K44,K33,Alpha2,[],[],true)
% plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% h = get(gcf,'children');
% for i = 1:length(h)
% 	if strcmp(h(i).Tag,'Colorbar')
% 		cbh = h(i);
% 		set(h(i),'YTick',-0.9:0.1:0.1)
% 	end
% end
savepcolor(gcf,'Figs/AlphaFB_needFFL')
xlabel('K_4')
ylabel('K_3')
title('\alpha_{FB}')


%
% FF/FB
%
figure
pcolor_contour(K44,K33,Alpha1,[],[],true)
% plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% h = get(gcf,'children');
% for i = 1:length(h)
% 	if strcmp(h(i).Tag,'Colorbar')
% 		cbh = h(i);
% 		set(h(i),'YTick',-0.9:0.1:0.1)
% 	end
% end
savepcolor(gcf,'Figs/AlphaFFFB_needFFL')
xlabel('K_4')
ylabel('K_3')
title('\alpha_{FFFB}')

%
% Now plot the difference b/w the two real parts. We will also add the
% stability contours.
%
dA = Alpha2 - Alpha1;
figure
pcolor_contour(K44,K33,dA,[-0.5:0.25:1],'r',true)
c = contourc(K44,K33,Alpha1,[0 0]);
[X1,Y1] = extractcontour(c);
plot(X1{1},Y1{1},'k');
c = contourc(K44,K33,Alpha2,[0 0]);
[X1,Y1] = extractcontour(c);
plot(X1{1},Y1{1},'k');
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
h = get(gcf,'children');
for i = 1:length(h)
	if strcmp(h(i).Tag,'Colorbar')
		cbh = h(i);
		set(h(i),'YTick',-0.6:0.2:1)
	end
end
savepcolor(gcf,'Figs/DeltaAlpha_needFFL')
xlabel('K_4')
ylabel('K_3')
title('\alpha_{FB}-\alpha_{FFFB}')



%}
