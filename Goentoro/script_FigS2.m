% script_FigS2
%
% This script is to plot the pcolor contour for Fig S2, which is Zpeak for
% the FF/FB model. We use some of the simulation results that we got when
% running a script embedded within script_Fig3.

clear
close all
load Mat/Fig3_all
% load Mat/Peaks_pointone_contours_FFFB X_realistic Y_realistic
load Mat/Peaks_pointone_contours XP YP Xzmax Yzmax
load Mat/C_pointone_contour XC YC

%
% Plotting pcolor and contour of Zpeak. This is for Fig. S2.
%
figure
pcolor_contour(K11,K22,Zpeak,[0.01 0.1 1],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(XP,YP,'Color',[1 0.5 0])
plot(XC,YC,'Color',[1 0.5 0])
plot(Xzmax,Yzmax,'w --')
% plot(X_realistic,Y_realistic,'w')
plot(K1_I_FBL,K2_I_FBL,'k','linewidth',1)
plot(K1_II_FBL,K2_II_FBL,'k','linewidth',1)
savepcolor(gcf,'Figs/Zmax_surface_FFFB')
xlabel('K1')
ylabel('K2')
title('Z_{max}')


% %
% % Get the zmax = 0.1 contour
% %
% h = get(gca,'children');
% for i = 1:length(h)
% 	if strcmp(h(i).DisplayName,'0.01')
% 		Xzmax = get(h(i),'Xdata')';
% 		Yzmax = get(h(i),'Ydata')';
% 		break
% 	end
% end
% nK2 = size(K12_NPA_FBL,1);
% K11_dense2 = logspace(-2,6,nK2)';
% K22_dense2 = logspace(-6,2,nK2)';
% Zpeak2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),Zpeak,...
% 	repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
% Vzmax = Zpeak2 >= 0.01;