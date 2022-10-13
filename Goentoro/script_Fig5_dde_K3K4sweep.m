% script_needFFL_dde_K48sweep
%
% This script is "goentoro" because we were using the goentoro model, with
% no delays. So, nominally, F0 = 1, then F1 = 10; So no need for
% constitutive synthesis.
%
% This script is "sweep" because we are sweeping K4 and K8, with K1 = 1;
%
% This script is designed to run code for a GRN that has both a negative
% feedback loop and a feedforward loop. The comparison here is between FB
% only and combined FF/FB.
%
% It is set up as follows.  The main output, Z, is repressed by W and Y.
% "Z" activates "W" . There is a disturbance, X, which activates "Z" and
% also "Y".  (This is the iffl: X -> Z; X -> Y -| Z).
%
% Each transcription/translation step is modeled as a delay differential
% equation, like in Longo et al., 2013.
%
% We will examine what "Z" does when the disturbance "X" undergoes a step
% increase from zero to one.  The negative feedback loop should kick in, so
% you should see an increase in "Z", followed by an increase in "W", then
% "Z" goes back down (to basal levels?).  If the feedback gain is too
% strong, or the deadtime is too long, we might see oscillation. How does
% this change if "Y" is knocked out vs present?
%
% When the iffl is present (so "Y" is not knocked out), then, under the
% same parameters, can we tune "Y" so that "Z" has the appropriate
% response?  So, eliminate or drastically dampen oscillations? OR, can we
% achieve a particular performance objective with "Y" present that is
% unavailable without it, no matter how we tune the neg fbk loop?

clear
close all
options = ddeset('RelTol',1e-6);
options2 = optimoptions('lsqcurvefit','display','off');

%
% Disturbance function
%
t0 = 0;
F0 = 1;
F1 = 10;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Param values
%
n = 2;
theta = [0.5 0.5 0.5];
K1 = 1; K2 = 0.1; 

y0 = F0.^n./(K1.^n + F0.^n);
y1 = F1.^n./(K1.^n + F1.^n);
RHS = y1'.^n - y0'.^n;
LHS = (1 + ((1./K2)*y0').^n)/F0.^n - (1 + ((1./K2)*y1').^n)/F1.^n;
% K3star = (repmat(RHS,nK2,1)./LHS).^(1/n);
K12 = (RHS./LHS).^(1/n);

tauz = 1; tauw = 1; tauy = 1;
tspan = [0 100];

nK = 50;
K33 = logspace(-3,1,nK);
K44 = logspace(-3,1,nK);

%% ========================================================================
% Running the set of simulations while varying K4,K8
% =========================================================================
% {

% figure
Zinitial = zeros(length(K33),length(K44),2);
Zpeak = zeros(length(K33),length(K44),2);
Tpeak = zeros(length(K33),length(K44),2);
Zfinal = zeros(length(K33),length(K44),2);
for i = 1:length(K33)
	for j = 1:length(K44)
		K3 = K33(i);
		K4 = K44(j);
		
		K5 = K1*K3; K6 = K2*K3; K7 = K12*K3;
		p = [tauz K1 K2 K12 K3 K4 tauw tauy];
		
		% =============================================================
		% FF/FB
		% =============================================================
		
		%
		% Initial and final conditions
		%
		[y0,z0,w0] = ftn_goentoro_ss(p,F0,n);
		[y1,z1,w1] = ftn_goentoro_ss(p,F1,n);
		
		
		%
		% FF/FB sim
		%
		ftnhand = @ftn_goentoro_dde;
		Y0 = [y0; z0; w0];
		soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
		t = soln.x';
		Y = soln.y';
		z = Y(:,2);
		Zinitial(i,j,1) = z0;
		[Zpeak(i,j,1),ipeak] = max(z);
		Tpeak(i,j,1) = t(ipeak) - t0;
		Zfinal(i,j,1) = z1;
		
% 		subplot(1,3,1)
% 		plot(t,z/z0)
% 		hold on
% 		subplot(1,3,2)
% 		plot(t,Y(:,3)/w0)
% 		hold on
% 		subplot(1,3,3)
% 		plot(t,Y(:,1)/y0)
% 		hold on
		
		
		% =============================================================
		% FB only
		% =============================================================

		%
		% Initial conditions
		%
		pFB = p; pFB(3) = Inf; pFB(4) = Inf;
		[y0FB,z0FB,w0FB] = ftn_goentoro_ss(pFB,F0,n);
		[y1FB,z1FB,w1FB] = ftn_goentoro_ss(pFB,F1,n);
		
		%
		% Running the FB only sim
		%
		Y0FB = [y0FB; z0FB; w0FB];
		solnfb = dde23(ftnhand,theta,Y0FB,tspan,options,theta,pFB,n,disthand);
		tFB = solnfb.x';
		YFB = solnfb.y';
		zFB = YFB(:,2);
		Zinitial(i,j,2) = z0FB;
		[Zpeak(i,j,2),ipeak] = max(zFB);
		Tpeak(i,j,2) = tFB(ipeak) - t0;
		Zfinal(i,j,2) = z1FB;
		
		
	end
	disp(['j = ',num2str(i),' out of ',num2str(nK)])
end


save Mat/2019-01-04_dde_K3K4sweep
%}

%%

%{
load Mat/2019-01-04_ss_K48sweep J_out
load Mat/2019-01-04_dde_K48sweep
load Mat/2018-11-30_eigval_of_J
close all

%
% FF/FB fold-change
%
F = abs(Zfinal-Zinitial)./Zinitial;
Alpha1 = interp1(J,Alpha,J_out(:,:,1));

% Stability contour for FF/FB
c = contourc(K44,K33,Alpha1,[0 0]);
[X1,Y1] = extractcontour(c);
K88star = X1{1}; K44star = Y1{1};
V1 = Alpha1 >= 0;

%
% FB only case
%
Alpha2 = interp1(J,Alpha,J_out(:,:,2));

% Stability contour for FB only
c = contourc(K44,K33,Alpha2,[0 0]);
[X1,Y1,C] = extractcontour(c);
K88star2 = X1{1}; K44star2 = Y1{1};
V2 = Alpha2 >= 0;

%
% Rel peak size, FB only
%
P = (Zpeak-Zinitial)./Zinitial;
DP = P(:,:,2); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
pcolor_contour(K44,K33,DP1,[0.8 0.9],'r',true)
shading interp
plot(K88star,K44star,'k','linewidth',2);
plot(K88star2,K44star2,'k','linewidth',2);
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/ratioP_needFFL')
xlabel('K_8')
ylabel('K_4')
title('P_{FB}')

%
% Zinitial, FF/FB
%
DP = Zinitial(:,:,1); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
pcolor_contour(K44,K33,DP1,[0.01 0.02],'r',true)
shading interp
plot(K88star,K44star,'k','linewidth',2);
plot(K88star2,K44star2,'k','linewidth',2);
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/ratioP_needFFL')
xlabel('K_8')
ylabel('K_4')
title('Z_0^{FFFB}')

%
% Abs peak size, FF/FB
%
DP = Zpeak(:,:,1); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
pcolor_contour(K44,K33,DP1,[0.02 0.03 0.04 0.05],'r',true)
shading interp
plot(K88star,K44star,'k','linewidth',2);
plot(K88star2,K44star2,'k','linewidth',2);
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/ratioP_needFFL')
xlabel('K_8')
ylabel('K_4')
title('Z_{max}^{FFFB}')


%
% Ratio of rel peak size, FF/FB to FF
%
P = (Zpeak-Zinitial)./Zinitial;
DP = P(:,:,1)./P(:,:,2); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
pcolor_contour(K44,K33,DP1,[1.4 1.5 1.6 1.7],'r',true)
shading interp
plot(K88star,K44star,'k','linewidth',2);
plot(K88star2,K44star2,'k','linewidth',2);
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/ratioP_needFFL')
xlabel('K_8')
ylabel('K_4')
title('P_{FFFB} / P_{FB}')

%
% Ratio of abs peak size, FF/FB to FF
%
DP = Zpeak(:,:,1)./Zpeak(:,:,2); 
DP1 = DP; DP11 = DP1(~(V1 | V2));
DP1(V1|V2) = min(DP11(:));
pcolor_contour(K44,K33,DP1,[0.1 0.2 0.3 0.4 0.5],'r',true)
shading interp
plot(K88star,K44star,'k','linewidth',2);
plot(K88star2,K44star2,'k','linewidth',2);
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/ratioP_needFFL')
xlabel('K_8')
ylabel('K_4')
title('Z_{max}^{FFFB} / Z_{max}^{FB}')

%}
