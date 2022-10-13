% script_needFFL_goentoro_dde_sweep
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

% clear
% close all

% {

%
% Changes to X
%
F0 = 1;
F1 = 10;

%
% Param values
%
n = 2;
K1 = 1; K2 = 0.1;

y0 = F0.^n./(K1.^n + F0.^n);
y1 = F1.^n./(K1.^n + F1.^n);
RHS = y1'.^n - y0'.^n;
LHS = (1 + ((1./K2)*y0').^n)/F0.^n - (1 + ((1./K2)*y1').^n)/F1.^n;
% K3star = (repmat(RHS,nK2,1)./LHS).^(1/n);
K12 = (RHS./LHS).^(1/n);

nK3 = 50;
nK4 = 50;
K33 = logspace(-3,1,nK3);
K44 = logspace(-3,1,nK4);

Zinitial = zeros(nK3,nK4,2);
Zfinal = zeros(nK3,nK4,2);
J_out = zeros(nK3,nK4,2);
for i = 1:nK3
	for j = 1:nK4
		K3 = K33(i);
		K4 = K44(j);
		
		p = [1 K1 K2 K12 K3 K4];
		
		% =============================================================
		% FF/FB
		% =============================================================
		
		%
		% Initial and final conditions
		%
		[y0,z0,w0] = ftn_goentoro_ss(p,F0,n);
		[y1,z1,w1] = ftn_goentoro_ss(p,F1,n);
		
		%
		% Recording
		%
		Zinitial(i,j,1) = z0;
		Zfinal(i,j,1) = z1;
		[~,J4,J5] = ftn_J4(K3,K1,K2,K3,K4,'K3',0,F1,n,true,[y1;z1;w1],K12);
		J_out(i,j,1) = J4*J5;
		
		
		% =============================================================
		% FB only
		% =============================================================

		%
		% Initial and final conditions
		%
		pFB = p; pFB(3) = Inf; pFB(4) = Inf;
		[y0FB,z0FB,w0FB] = ftn_goentoro_ss(pFB,F0,n);
		[y1FB,z1FB,w1FB] = ftn_goentoro_ss(pFB,F1,n);
		
		%
		% Recording
		%
		Zinitial(i,j,2) = z0FB;
		Zfinal(i,j,2) = z1FB;
		[g,J4,J5] = ftn_J4(K3,K1,0,K3,K4,'K3',0,F1,n,false,[0;z1FB;w1FB]);
		J_out(i,j,2) = J4*J5;
		
	end
end

save Mat/2019-01-04_ss_K3K4sweep J_out K33 K44 Zinitial Zfinal K1 K2 K12 n F0 F1
%}

%% ========================================================================
% Plotting.
% =========================================================================

%{
% clear
load Mat/2019-01-04_ss_K3K4sweep J_out K33 K44
load Mat/2018-11-30_eigval_of_J J Alpha

close all

%
% Now that we have J_out, we can convert those numbers into eigenvalues.
%
Alpha1 = interp1(J,Alpha,J_out(:,:,1));
Alpha2 = interp1(J,Alpha,J_out(:,:,2));

pcolor_contour(K44,K33,Alpha1,[],[],true)
% pcolor_contour(K44,K33,Alpha1,[0.1 0 -0.5 -0.9],{'r','k','r','r'},true)
xlabel('K_4')
ylabel('K_3')
title('\alpha_{FFFB}')
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)

pcolor_contour(K44,K33,Alpha2,[],[],true)
% pcolor_contour(K44,K33,Alpha2,[0.1 0 -0.5 -0.9],{'r','k','r','r'},true)
xlabel('K_4')
ylabel('K_3')
title('\alpha_{FB}')
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)

%
% Now plot the difference b/w the two real parts. We will also add the
% stability contours.
%
pcolor_contour(K44,K33,Alpha2-Alpha1,[],[],true)
% pcolor_contour(K44,K33,Alpha2-Alpha1,[0.5 0 -0.5],{'r','k','r'},true)
xlabel('K_4')
ylabel('K_3')
title('\alpha_{FB}-\alpha_{FFFB}')
c = contourc(K44,K33,Alpha1,[0 0]);
[X1,Y1] = extractcontour(c);
plot(X1{1},Y1{1},'k');
c = contourc(K44,K33,Alpha2,[0 0]);
[X1,Y1,C] = extractcontour(c);
plot(X1{1},Y1{1},'k');
plot([0.1 0.1 0.5 0.5 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)


%}















