% script_FigS5
%
% This script is designed to find the DeltaK12NPAFFFB/DeltaK12NPAFF for
% different values of epsilon.

clear
close all
options = ddeset('RelTol',1e-6);

%
% Param values
%
x0 = 1;
x1 = 10;
n = 1;
vareps = 0.05;
Vareps = [0.01 0.1];
epsname = {'pointohone','pointone'};

thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
tauz = 1; tauw = 1; tauy = 1;
K3 = 0.1; K4 = 0.1;

%
% Disturbance function
%
t0 = 0;
tspan = [t0 100];
disthand = @(t)ftn_unitstep(t,t0,x0,x1);




%% ========================================================================
% First, we will calculate the Delta(K12)/K12 for the situation with
% FF only. Then the sweep for the one with neg FBK.
% =========================================================================
% {

for iDU = 1:length(Vareps)
	vareps = Vareps(iDU);

nK = 500;
K11 = logspace(-2,6,nK)'; K111 = repmat(K11',nK,1);
K22 = logspace(-6,2,nK)';

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
RHS = K11'.*(y1' - y0'/(1 + vareps));
LHS = 1/(1 + vareps)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAplus = repmat(RHS,nK,1)./LHS;

% K12_NPAminus
vareps1 = -vareps;
RHS = K11'.*(y1' - y0'/(1 + vareps1));
LHS = 1/(1 + vareps1)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAminus = repmat(RHS,nK,1)./LHS;

K1lower_minus = vareps1/(1/x1 - 1/x0*(1 + vareps1));
K12_NPAminus(K11 <= K1lower_minus) = 0;
K12_NPAminus(K12_NPAminus < 0) = 0;

Y = (K12_NPAplus-K12_PA)./K12_PA;
Yminus = (K12_PA-K12_NPAminus)./K12_PA;
DeltaK = Y + Yminus;

%
% Regions I and II
%
script_regionsI_and_II


%
% creating K12_NPAplus and minus arrays
%
script_Fig3_FBL_K12NPA_K1K2sweep
Y = (real(K12_NPA_FBL(:,:,1)) - K12_PA)./K12_PA; 
Yminus = (K12_PA - real(K12_NPA_FBL(:,:,2)))./K12_PA;
DeltaK_FBL = Y + Yminus;


%% ========================================================================
% Now plot
% =========================================================================
% {
load Mat/Peaks_pointone_contours_FFFB X_realistic Y_realistic

% %
% % Just for FFL only
% %
% figure
% pcolor_contour(K11,K22,DeltaK,[0.15 0.3 0.67 5],'r',true,0,5)%,0.1223,5)
% % set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
% plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
% plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
% plot(X_realistic,Y_realistic,'w')
% % h = get(gcf,'children');
% % for i = 1:length(h)
% % 	if strcmp(h(i).Tag,'Colorbar')
% % 		cbh = h(i);
% % 		set(h(i),'YTick',1:0.25:2)
% % 	end
% % end
% % savepcolor(gcf,'Figs/Ratio_of_DeltaK12_difteps')
% xlabel('K1')
% ylabel('K2')
% title('\Delta K_{12}^{NPA}/K_{12}^{NPA}')


RDK = DeltaK_FBL./DeltaK; % ratio of DeltaK's
% RDK = DeltaK; % ratio of DeltaK's
figure
pcolor_contour(K11,K22,RDK,[1 1.1 1.25 1.5 1.75 2],'r',true,0.99,2)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
plot(X_realistic,Y_realistic,'w')
h = get(gcf,'children');
for i = 1:length(h)
	if strcmp(h(i).Tag,'Colorbar')
		cbh = h(i);
		set(h(i),'YTick',1:0.25:2)
	end
end
savepcolor(gcf,['Figs/Ratio_of_DeltaK12_eps_',epsname{iDU}])
xlabel('K1')
ylabel('K2')
title('\Delta K_3^{NPA,FBL}/\Delta K_3^{NPA}')

end

%}

