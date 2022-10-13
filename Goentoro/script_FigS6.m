% script_FigS6
%
% This script, for Fig. S6, is to find the parameter regimes that give
% different phenotypes. All of these calculations are for the FF only model
% with PA (so K12 is finely tuned).
%
% We will look at three phenotypes: (1) pulse generation, (2) response
% acceleration, and (3) fold change detection.
%
% We will not look at perfect adaptation, since that has already been
% dealty with, obviously, since that's the main focus of the paper.
%
% For pulse generation, that requires both P > 0.1 AND Q > 0.1. The P > 0.1
% contour has already been found. But we don't have to do Q, because we are
% running all of these at PA, so P = Q (cause Zinit = Zfinal).
%
% For response acceleration, we will
% 
% For FCD, we will use Goentoro et al's formulation, and say FCD happens
% when the two peaks (one for 1 -> 10, the other for 5 -> 50) are within
% 10% of each other.
%

clear
close all
options = ddeset('RelTol',1e-6);
load Mat/Fig2_peak_K1K2sweep
load Mat/C_pointone_contour XC YC VC1
load Mat/Peaks_pointone_contours XP YP Xzmax Yzmax VP Vzmax

K11_dense = logspace(log10(K11(1)),log10(K11(end)),size(VC1,1))';
K22_dense = logspace(log10(K22(1)),log10(K22(end)),size(VC1,1))';

%
% Param values
%
x0 = 1;
x1 = 10;
n = 1;
vareps = 0.1;

thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
tauz = 1; tauw = 1; tauy = 1;


%
% Disturbance function
%
t0 = 0;
tspan = [t0 20];
disthand = @(t)ftn_unitstep(t,t0,x0,x1);


%
% Create realistic region
%
V_realistic = VP & VC1 & Vzmax;
[~,J_realistic,I_realistic] = traceobject(V_realistic);
X_realistic = K11_dense(J_realistic);
Y_realistic = K22_dense(I_realistic);





%% ========================================================================
% Finally, we will calculate the peak for a FC of 10 when going from 5 to
% 50.
% =========================================================================
% {

x0_2 = 5;
x1_2 = 50;

y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
RHS = y1' - y0';
LHS = (1 + (1./K22)*y0')/x0 - (1 + (1./K22)*y1')/x1;
K12_PA = repmat(RHS,nK,1)./LHS;

disthand2 = @(t)ftn_unitstep(t,t0,x0_2,x1_2);
% script_FigS6_FFL_dde_K1K2sweep
load Mat/FigS6_peak_K1K2sweep
% load Mat/Zfinal_21
% Z_final2 = Zfinal_21;


O = (Zpeak_2 - Zpeak)./Zpeak;
% P = (Zpeak - Zinitial)./Zinitial;
% P_2 = (Zpeak_2 - Zinitial_2)./Zinitial_2;
% I = (Zinitial_2-Zinitial)./Zinitial;
% F = (Zfinal_2-Zfinal)./Zfinal;
% RP = P./P_2;

figure
pcolor_contour(K11,K22,smooth2b(abs(O),5,5),[0.05 0.1 0.2 0.3 0.4],'r',true)%,0.99,2)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]),'fontsize',24)
plot(X_realistic,Y_realistic,'w')
% h = get(gcf,'children');
% for i = 1:length(h)
% 	if strcmp(h(i).Tag,'Colorbar')
% 		cbh = h(i);
% 		set(h(i),'YTick',1:0.25:2)
% 	end
% end
savepcolor(gcf,'Figs/FCD')
xlabel('K1')
ylabel('K2')
title('I')


%}


%% ========================================================================
% We have some simulations that are not quite working out...for some
% reason, our simulations predict that Region II is the worst place for
% FCD, even though that's supposed to be the best place. I am wondering if
% that's because of our choice of K12?
% =========================================================================
%{

%
% Param values: run 1
%
n = 1;
x0 = 1;
x1 = 10;
theta = [0.5 0.5 0.5];
% K1 = 1060; K2 = 1.012e-5; K3 = Inf; K4 = Inf;
K1 = K11(126); K2 = K22(26); K3 = Inf; K4 = Inf;
tauz = 1; tauw = 1; tauy = 1;

% In case you're looking for PA
y0 = x0./(K1 + x0);
y1 = x1./(K1 + x1);
RHS = y1' - y0';
LHS = (1/x0-1/x1) - (1./K2)*(1./(K1' + x1) - 1./(K1' + x0));
K12 = RHS./LHS;


p = [tauz K1 K2 K12 K3 K4 tauw tauy];
yesplot = true;
figure
ftn_rundde(p,theta,x0,x1,n,yesplot,[0 20]);


%
% Param values: run 2
%
n = 1;
x0 = 5;
x1 = 50;
theta = [0.5 0.5 0.5];
K1 = K11(126); K2 = K22(26); K3 = Inf; K4 = Inf;
tauz = 1; tauw = 1; tauy = 1;

% In case you're looking for PA
y0 = x0./(K1 + x0);
y1 = x1./(K1 + x1);
RHS = y1' - y0';
LHS = (1/x0-1/x1) - (1./K2)*(1./(K1' + x1) - 1./(K1' + x0));
% K12 = RHS./LHS;


p = [tauz K1 K2 K12 K3 K4 tauw tauy];
yesplot = true;
hold on
ftn_rundde(p,theta,x0,x1,n,yesplot,[0 10]);





%}


%% ========================================================================
% First, we will calculate Q. Oops, Q = P because of the PA constraint.
% =========================================================================
%{



Q = (Zpeak - Zfinal)./Zfinal;
figure
pcolor_contour(K11,K22,Q,[0.1 0.5 1 3],'r',true)%,0.99,2)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(X_realistic,Y_realistic,'w')
% h = get(gcf,'children');
% for i = 1:length(h)
% 	if strcmp(h(i).Tag,'Colorbar')
% 		cbh = h(i);
% 		set(h(i),'YTick',1:0.25:2)
% 	end
% end
% savepcolor(gcf,['Figs/Ratio_of_DeltaK12_eps_',epsname{iDU}])
xlabel('K1')
ylabel('K2')
title('Q')


%}


