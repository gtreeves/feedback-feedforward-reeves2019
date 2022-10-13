% script_allconstraints
%
% This script is designed to find the region of K1,K2, K4,K8 param space
% that satisfies all important constraints.
% (1) PA: this is satisfied by choosing the correct K3.
% (2) z0 > 0.01 (realistic constraint)
% (3) zmax > 0.05 (need to be guaranteed of a minimum abs height peak)
% (4) P > 0.1 (need good signal; NOTE: SAR is undefined)
% (5) stable (of course)
% (6) minimum settling time = 50?

clear
close all
options = ddeset('RelTol',1e-6);

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
K1 = 1; K2 = 0.3; K3 = 0.3; K4 = 0.01; K8 = 0.01;
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 100];
% L = 1000;
% tstep = linspace(tspan(1),tspan(2),L)';
% [~,i50] = min(abs(tstep - tspan(2)));

nK4 = 50;
nK8 = 50;
K22 = [0.01 0.1 1]';
nK2 = length(K22);
K44 = logspace(-2,0,nK4);
K88 = logspace(-2,0,nK8);
% K88 = 0.1;

K11 = 1; nK1 = length(K11);
y0 = F0.^n./(K11.^n + F0.^n);
y1 = F1.^n./(K11.^n + F1.^n);
RHS = y1'.^n - y0'.^n;
% RHS = (F1./(K11' + F1)).^n - (F0./(K11' + F0)).^n;
% LHS = (1/F0^n-1/F1^n) - (1./K22)*(1./(K11' + F1) - 1./(K11' + F0));
LHS = (1 + ((1./K22)*y0').^n)/F0.^n - (1 + ((1./K22)*y1').^n)/F1.^n;
K3star = (repmat(RHS,nK2,1)./LHS).^(1/n);

K5 = K1*K4; K6 = K2*K4; K7 = K3*K4;


%% ========================================================================
% Now we will run a dde sweep of K4,K8, embedded inside of a sweep of K2.
% =========================================================================
% {

load Mat/2018-11-30_eigval_of_J

Zinitial = NaN(nK4,nK8,nK2);
Zpeak = NaN(nK4,nK8,nK2);
Tsettle = NaN(nK4,nK8,nK2);
Alpha1 = NaN(nK4,nK8,nK2);
ftnhand = @ftn_goentoro_dde;

for ii = 1:nK2
	
	K2 = K22(ii);
	jj = 1;
	K1 = K11(jj);
	K3 = K3star(ii,jj);
	
	for i = 1:nK4
		for j = 1:nK8
			K4 = K44(i);
			K8 = K88(j);
			
			K5 = K1*K4; K6 = K2*K4; K7 = K3*K4;
			p = [tauz K1 K2 K3 K4 K8 tauw tauy];
			
			% =============================================================
			% FF/FB
			% =============================================================
			
			%
			% Initial and final conditions
			%
			[y0,z0,w0] = ftn_goentoro_ss(p,F0,n);
			[y1,z1,w1] = ftn_goentoro_ss(p,F1,n);
			
			
			%
			% Stability
			%
			% This function, "ftn_J4", calculates our two partials J4 and
			% J5 given values of K1,K2,K3,K4,K8. (If no K3 is given, then c
			% = 1).
			%
			[~,J4,J5] = ftn_J4(K4,K1,K2,K4,K8,'K4',0,F1,n,true,[y1;z1;w1],K3);
			J_out = J4*J5;
			Alpha1(i,j,ii) = interp1(J,Alpha,J_out); % if neg, then stable
			
			if Alpha1(i,j,ii) < 0
			
				%
				% Combined FF/FB sim
				%
				Y0 = [y0; z0; w0];
				soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
				t = soln.x';
				Y = soln.y';
				z = Y(:,2);
				Zinitial(i,j,ii) = z0;
				[Zpeak(i,j,ii),ipeak] = max(z);
				
				k = find(abs(z-z0)/z0 > 0.05);
				k(k < ipeak) = []; k = k(1);
				Tsettle(i,j,ii) = t(k) - t0;
			
			end
			
			
		end
	end
	disp(['ii = ',num2str(ii),' out of ',num2str(nK2)])
end

save Mat/2018-12-25_allconstraints

%}

%% ========================================================================
% Now that we have run a dde sweep of K4,K8, embedded inside of a sweep of
% K2, we can plot the pcolor contours.
% =========================================================================
% {

close all
load Mat/2018-12-25_allconstraints

P = (Zpeak - Zinitial)./Zinitial;
for ii = 1:nK2
	
	pcolor_contour(K88,K44,Alpha1(:,:,ii),Alpha1(:,:,ii),0,'r',true)
	xlabel('K_8')
	ylabel('K_4')
	title(['\alpha , K_2 = ',num2str(K22(ii))])
	
	pcolor_contour(K88,K44,Zinitial(:,:,ii),Zinitial(:,:,ii),0,'r',true)
	xlabel('K_8')
	ylabel('K_4')
	title(['z_0, K_2 = ',num2str(K22(ii))])
	
	pcolor_contour(K88,K44,Zpeak(:,:,ii),Zpeak(:,:,ii),0,'r',true)
	xlabel('K_8')
	ylabel('K_4')
	title(['z_{max}, K_2 = ',num2str(K22(ii))])
	
	pcolor_contour(K88,K44,P(:,:,ii),P(:,:,ii),0,'r',true)
	xlabel('K_8')
	ylabel('K_4')
	title(['P, K_2 = ',num2str(K22(ii))])
	
	pcolor_contour(K88,K44,Tsettle(:,:,ii),Tsettle(:,:,ii),0,'r',true)
	xlabel('K_8')
	ylabel('K_4')
	title(['t_{settle}, K_2 = ',num2str(K22(ii))])
	
end











