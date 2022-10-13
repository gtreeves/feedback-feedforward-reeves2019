% script_Fig5
%
% This script is designed to run all simulations, etc, to get Fig 5 of
% Reeves 2019 (FF/FB).
%
% What do we want in Fig 5? This fig is for a FB system as baseline, and
% how it could benefit from adding FF control. 
%
% Part (a): Difference between 
%
% Maybe hold the Adaptation, Peak, Other peak, and SNR pcolor contours for
% Fig 5?

clear
close all
options = ddeset('RelTol',1e-6);


%% ========================================================================
% No we'll do ratio of F, of P, of O, of SAR, and SAR_alt.
% =========================================================================
% {


load Mat/2019-01-04_ss_K3K4sweep J_out
load Mat/2019-01-04_dde_K3K4sweep
load Mat/2018-11-30_eigval_of_J
close all

%
% FF/FB fold-change
%
F = abs(Zfinal-Zinitial)./Zinitial;
Alpha1 = interp1(J,Alpha,J_out(:,:,1));
Beta1 = interp1(J,Beta,J_out(:,:,1));

% Stability contour for FF/FB
c = contourc(K44,K33,Alpha1,[0 0]);
[X1,Y1] = extractcontour(c);
K44star = X1{1}; K33star = Y1{1};
V1 = Alpha1 >= 0;

%
% FB only case
%
Alpha2 = interp1(J,Alpha,J_out(:,:,2));
Beta2 = interp1(J,Beta,J_out(:,:,2));

% Stability contour for FB only
c = contourc(K44,K33,Alpha2,[0 0]);
[X1,Y1,C] = extractcontour(c);
K44star2 = X1{1}; K33star2 = Y1{1};
V2 = Alpha2 >= 0;

%
% Rel peak size, FB only
%
P = (Zpeak-Zinitial)./Zinitial;
DP = P(:,:,2); 
DP1 = DP; DP11 = DP1(~(V2));
DP1(V2) = min(DP11(:));
figure
pcolor_contour(K44,K33,{DP1 DP},[],'r',true,[],[])%,0.7:0.05:0.95)
shading interp
plot(K44star,K33star,'k','linewidth',2);
plot(K44star2,K33star2,'k','linewidth',2);
% plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
savepcolor(gcf,'Figs/P_FBonly_needFFL')
xlabel('K_4')
ylabel('K_3')
title('P_{FB}')


%
% Ratio of peak size, FF/FB to FF
%
P = (Zpeak-Zinitial)./Zinitial;
DP = P(:,:,1)./P(:,:,2); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
figure
pcolor_contour(K44,K33,{DP1 DP},[],'r',true,[],[])%,1.3:0.1:1.8)
shading interp
plot(K44star,K33star,'k','linewidth',2);
plot(K44star2,K33star2,'k','linewidth',2);
% plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% h = get(gcf,'children');
% for i = 1:length(h)
% 	if strcmp(h(i).Tag,'Colorbar')
% 		cbh = h(i);
% 		set(h(i),'YTick',0.48:0.02:0.64)
% 	end
% end
savepcolor(gcf,'Figs/ratioP_needFFL')
xlabel('K_4')
ylabel('K_3')
title('P_{FFFB}/P_{FB}')

%
% Abs peak size, FF/FB
%
DP = Zpeak(:,:,1); 
DP1 = DP; DP11 = DP1(~(V1));
DP1(V1) = min(DP11(:));
figure
pcolor_contour(K44,K33,{DP1 DP},[],'r',true,[],[])%,0.02:0.01:0.06)
shading interp
plot(K44star,K33star,'k','linewidth',2);
plot(K44star2,K33star2,'k','linewidth',2);
% plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
savepcolor(gcf,'Figs/zmaxFFFB_needFFL')
xlabel('K_4')
ylabel('K_3')
title('Z_{max}^{FFFB}')


%
% Ratio of peak size, FF/FB
%
DP = Zpeak(:,:,1)./Zpeak(:,:,2); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
figure
pcolor_contour(K44,K33,{DP1 DP},[],'r',true,[],[])%,0.1:0.1:0.6)
shading interp
plot(K44star,K33star,'k','linewidth',2);
plot(K44star2,K33star2,'k','linewidth',2);
% plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
savepcolor(gcf,'Figs/ratioZmax_needFFL')
xlabel('K_4')
ylabel('K_3')
title('Z_{max}^{FFFB}/Z_{max}^{FB}')




%}