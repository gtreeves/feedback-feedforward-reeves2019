% script_Fig2
%
% This script is designed to run all simulations, etc, to get Fig 2 of
% Reeves 2019 (FF/FB).

clear
close all
options = ddeset('RelTol',1e-6);

vareps = 0.05;

%% ========================================================================
% First, we find the value of K12 for which f = eps = 0.05. Then for f =
% -0.05. We will do this for a particular choice of K1 and K2. This is so
% we can run a simulation and show what it looks like.
% =========================================================================
% {

K1 = 1; K2 = 0.1; 

%
% K12_PA
%
x0 = 1;
x1 = 10;
y0 = x0./(K1 + x0);
y1 = x1./(K1 + x1);
RHS = y1 - y0;
LHS = (1/x0-1/x1) - (1./K2)*(1./(K1 + x1) - 1./(K1 + x0));
K12 = RHS/LHS;



%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,x0,x1);

%
% Param values
%
n = 1;
thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 20];
K3 = 1; K4 = Inf;
p = [tauz K1 K2 K12 K3 K4 tauw tauy];

%
% Initial conditions
%
y0 = x0./(K1 + x0);
z0 = x0/K1./(1 + x0/K1 + y0/K2 + x0*y0/K12);
w0 = 0;

%
% Perfect adaptation FF only sim
%
ftnhand = @ftn_goentoro_dde;
Y0 = [y0; z0; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
t = soln.x';
Y = soln.y';
z = Y(:,2);

%
% creating K12_NPAplus
%
RHS = K1.*(y1 - y0/(1 + vareps));
LHS = 1/(1 + vareps)*(1 + K1/x0 + (1./K2/x0)*(K1.*y0)) - ...
	(1 + K1/x1 + (1./K2/x1)*(K1.*y1));
K12_NPAplus = RHS/LHS;
p(4) = K12_NPAplus;

%
% Near perfect adaptation (plus), FF only sim
%
z0plus = x0/K1./(1 + x0/K1 + y0/K2 + x0*y0/K12_NPAplus);
Y0(2) = z0plus;
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
tplus = soln.x';
Y = soln.y';
zplus = Y(:,2);


%
% creating K12_NPAminus
%
vareps1 = -vareps;
RHS = K1.*(y1 - y0/(1 + vareps1));
LHS = 1/(1 + vareps1)*(1 + K1/x0 + (1./K2/x0)*(K1.*y0)) - ...
	(1 + K1/x1 + (1./K2/x1)*(K1.*y1));
K12_NPAminus = RHS/LHS;
p(4) = K12_NPAminus;

%
% Near perfect adaptation (plus), FF only sim
%
z0minus = x0/K1./(1 + x0/K1 + y0/K2 + x0*y0/K12_NPAminus);
Y0(2) = z0minus;
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
tminus = soln.x';
Y = soln.y';
zminus = Y(:,2);


%
% Plotting
%
figure
plot(t,z,tplus,zplus,tminus,zminus,'linewidth',2)
set(gca,'fontsize',24)%,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
print(gcf,'Figs/NPA_simulation_FFonly.eps','-depsc')
print(gcf,'Figs/NPA_simulation_FFonly.jpg','-djpeg','-r150')

figure
plot(t,z/z0,tplus,zplus/z0plus,tminus,zminus/z0minus,'linewidth',2)
hold on
plot(xlim,[1 1]*1.05,'k:',xlim,[1 1]*0.95,'k:')
set(gca,'fontsize',24)%,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
print(gcf,'Figs/NPA_simulation_FFonly_norm.eps','-depsc')
print(gcf,'Figs/NPA_simulation_FFonly_norm.jpg','-djpeg','-r150')



%}





%% ========================================================================
% Next, we plot the surface of K12 for PA. We will do this for varying K1
% and K2 all over the map. 
% =========================================================================
% {
nK = 200;
K11 = logspace(-2,6,nK)';
K22 = logspace(-6,2,nK)';

%
% K3_PA
%
x0 = 1;
x1 = 10;
y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
% LHS = (1/F0-1/F1) - (1./(K1 + F1) - 1./(K1 + F0))*(1./K2');
RHS = y1' - y0';
LHS = (1/x0-1/x1) - (1./K22)*(1./(K11' + x1) - 1./(K11' + x0));
K12_PA = repmat(RHS,nK,1)./LHS;

%
% Boundaries of Regions I and II
%
RHS = 1/(1 + vareps)*(K11.*y0/x0) - (K11.*y1/x1);
LHS = (1 + K11/x1) - 1/(1 + vareps)*(1 + K11/x0);
K2_I = RHS./LHS;

v = K2_I > 0;
[V,v_val] = repeatcheck(v);
k = find(v_val == 1);

K2_II = K2_I(V{k(2)});
K2_I = K2_I(V{k(1)});

K1_I = K11(V{k(1)});
K1_II = K11(V{k(2)});
K1upper = vareps/(1/x0 - (1 + vareps)/x1);
K1lower = (x1 - (1 + vareps)*x0)/vareps;

% %
% % pcolor (w/contours) of K12^{PA} vs K1,K2
% %
% pcolor_contour(K11,K22,log10(K12_PA),[-7 -6 -5 -4 -3 -2 -1 log10(0.5)],'r',true)
% set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
% plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
% plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
% 
% savepcolor(gcf,'Figs/K12_PA_surface')
% xlabel('K1')
% ylabel('K2')
% title('K_3^{PA}')


%
% pcolor (w/contours) of C = K3^{PA}/(K1*K2)
%
C = K12_PA./repmat(K11',nK,1)./repmat(K22,1,nK);
figure('paperpositionmode','auto')
pcolor_contour(K11,K22,log10(C),[-12:-1 log10(0.95)],'r',true,-6,0)%
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
% plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
savepcolor(gcf,'Figs/cooperativity_surface')
xlabel('K1')
ylabel('K2')
title('C = K_{12}^{PA}/(K_1K_2)')

%
% Get the C = 0.01 contour
%
h = get(gca,'children');
for i = 1:length(h)
	if strcmp(h(i).DisplayName,'-2')
		XC = get(h(i),'Xdata')';
		YC = get(h(i),'Ydata')';
		break
	end
end

nK1 = 1000;
K11_dense1 = logspace(-2,6,nK1)';
K22_dense1 = logspace(-6,2,nK1)';
nK2 = 500;
K11_dense2 = logspace(-2,6,nK2)';
K22_dense2 = logspace(-6,2,nK2)';
C1 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),C,...
	repmat(K11_dense1',nK1,1),repmat(K22_dense1,1,nK1));
VC1 = C1 >= 0.01;
C2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),C,...
	repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
VC = C2 >= 0.01;

save Mat/C_pointone_contour XC YC VC VC1
%}


%% ========================================================================
% Now we run a parameter sweep to get the peak value parameter P while
% varying values of K1,K2, at K3_PA, for FFL only. It turns out that this
% parameter is very weak for Region I. So that's a no-go region.
% =========================================================================
% {
nK = 200;
K11 = logspace(-2,6,nK)';
K22 = logspace(-6,2,nK)';
x0 = 1;
x1 = 10;

%
% K3_PA
%
y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
y00 = repmat(y0',nK,1);
y11 = repmat(y1',nK,1);
% RHS = y11 - y00;
% LHS = (1 + (1./K22)*y0')/F0 - (1 + (1./K22)*y1')/F1;
% K3_PA = RHS./LHS;

%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,x0,x1);

%
% Param values
%
n = 1;
thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
K1 = 1; K2 = 0.1; K12 = 0.3; K3 = 0.01; K4 = 0.01;
K5 = K1*K3; K6 = K2*K3; K7 = K12*K3;
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 100];

%
% All (slow) simulation commands in the following script, which is
% commented out now, since the sims were already run and stored in the Mat
% file below.
%
% script_Fig2_FFL_dde_K1K2sweep
load Mat/Fig2_peak_K1K2sweep

%
% Boundary of Region I (and II, too)
%
% vareps = 0.05;
nK2 = 1000;
K111 = logspace(-2,6,nK2)';
y0 = x0./(K111 + x0);
y1 = x1./(K111 + x1);
RHS = 1/(1 + vareps)*(K111.*y0/x0) - (K111.*y1/x1);
LHS = (1 + K111/x1) - 1/(1 + vareps)*(1 + K111/x0);
K2_I = RHS./LHS;
v = K2_I > 0;
[V,v_val] = repeatcheck(v);
k = find(v_val == 1);

K2_II = K2_I(V{k(2)});
K2_I = K2_I(V{k(1)});

K1_I = K111(V{k(1)});
K1_II = K111(V{k(2)});
K1upper = vareps/(1/x0 - (1 + vareps)/x1);
K1lower = (x1 - (1 + vareps)*x0)/vareps;


%
% Plotting pcolor and contour of P.
%
P = (Zpeak - Zinitial)./Zinitial;
% P1 = smooth2a(P,3,3);
% P_RH = P(:,round(nK/2):end);
% P_RH1 = smooth2a(P_RH,15,15);
% P1 = [P(:,1:round(nK/2)-1) P_RH1];

figure('paperpositionmode','auto')
pcolor_contour(K11,K22,P,[0.01 0.1 1 1.95 2.05 3],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
savepcolor(gcf,'Figs/Peak_surface')
xlabel('K1')
ylabel('K2')
title('P = (Z_{max}-Z_0)/Z_0')

%
% Get the P = 0.1 contour
%
h = get(gca,'children');
for i = 1:length(h)
	if strcmp(h(i).DisplayName,'0.1')
		XP = get(h(i),'Xdata')';
		YP = get(h(i),'Ydata')';
		break
	end	
end
nK2 = 1000;
K11_dense2 = logspace(-2,6,nK2)';
K22_dense2 = logspace(-6,2,nK2)';
P2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),P,...
	repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
VP = P2 >= 0.1;

%
% Plotting pcolor and contour of Zpeak.
%
figure('paperpositionmode','auto')
pcolor_contour(K11,K22,Zpeak,[0.0001 0.001 0.01 0.1 1],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(XP,YP,'color',[1 0.5 0])
plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
savepcolor(gcf,'Figs/Zmax_surface')
xlabel('K1')
ylabel('K2')
title('Z_{max}')

%
% Get the zmax = 0.1 contour
%
h = get(gca,'children');
for i = 1:length(h)
	if strcmp(h(i).DisplayName,'0.01')
		Xzmax = get(h(i),'Xdata')';
		Yzmax = get(h(i),'Ydata')';
		break
	end
end

nK2 = 1000;
K11_dense2 = logspace(-2,6,nK2)';
K22_dense2 = logspace(-6,2,nK2)';
Zpeak2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),Zpeak,...
	repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
Vzmax = Zpeak2 >= 0.01;

save Mat/Peaks_pointone_contours XP YP Xzmax Yzmax VP Vzmax

%}


%% ========================================================================
% Next, we find the surface of K12 for which f = eps = 0.05. Then for f =
% -0.05. We will do this for varying K1 and K2 all over the map.
% =========================================================================
% {
nK = 1000;
K11 = logspace(-2,6,nK)'; K111 = repmat(K11',nK,1);
K22 = logspace(-6,2,nK)';

load Mat/C_pointone_contour XC YC VC1
load Mat/Peaks_pointone_contours XP YP Xzmax Yzmax VP Vzmax

%
% K12_PA
%
x0 = 1;
x1 = 10;
y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
RHS = y1' - y0';
LHS = (1/x0-1/x1) - (1./K22)*(1./(K11' + x1) - 1./(K11' + x0));
K12_PA = repmat(RHS,nK,1)./LHS;

%
% creating K12_NPAplus and minus arrays
%
% vareps = 0.05;
RHS = K11'.*(y1' - y0'/(1 + vareps));
LHS = 1/(1 + vareps)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAplus = repmat(RHS,nK,1)./LHS;

% Boundary of Region I (and II, too)
RHS = 1/(1 + vareps)*(K11.*y0/x0) - (K11.*y1/x1);
LHS = (1 + K11/x1) - 1/(1 + vareps)*(1 + K11/x0);
K2_I = RHS./LHS;

v = K2_I > 0;
[V,v_val] = repeatcheck(v);
k = find(v_val == 1);

K2_II = K2_I(V{k(2)});
K2_I = K2_I(V{k(1)});

K1_I = K11(V{k(1)});
K1_II = K11(V{k(2)});
K1upper = vareps/(1/x0 - (1 + vareps)/x1);
K1lower = (x1 - (1 + vareps)*x0)/vareps;

% K12_NPAminus
vareps1 = -vareps;
RHS = K11'.*(y1' - y0'/(1 + vareps1));
LHS = 1/(1 + vareps1)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAminus = repmat(RHS,nK,1)./LHS;

K1lower_minus = vareps1/(1/x1 - 1/x0*(1 + vareps1));
K12_NPAminus(K11 <= K1lower_minus) = 0;
K12_NPAminus(K12_NPAminus < 0) = 0;


%
% Plotting
%
Y = (K12_NPAplus-K12_PA)./K12_PA;
Yminus = (K12_PA-K12_NPAminus)./K12_PA;
DeltaK = Y + Yminus;
% K12_NPAminus1 = K12_NPAminus;
% K12_NPAminus1(K12_NPAminus < 0) = 0;
figure('paperpositionmode','auto')
pcolor_contour(K11,K22,DeltaK,[0.15 0.3 0.67 5],'r',true,0,5)%,0.1223,5)
% pcolor_contour(K11,K22,K12_NPAminus1,[0 0.15 0.3 1 5 10],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
plot(XC,YC,'w',XP,YP,'w',Xzmax,Yzmax,'w')
savepcolor(gcf,'Figs/DeltaK12_NPA_FFLonly')
xlabel('K1')
ylabel('K2')
title('\Delta K_{12}^{NPA}/K_{12}^{PA}')

K11_dense = K11;
K22_dense = K22; % these are denser versions, something we can afford in
% this cell, since the calculation is so fast.
save Mat/DeltaK_FFLonly_dense K11_dense K22_dense DeltaK

%
% Create mathematical definition of realistic region
%
V = VP & VC1 & Vzmax;
min(DeltaK(V))
max(DeltaK(V))

%}




















