% script_goentoro
%
% This script initiates the equations from Goentoro et al., 2009 to
% replicate some of her results. In particular, we will investigate the
% region of the parameter space in which you get perfect or near-perfect
% adaptation.
%
% To have the fold-change detection, she assumed that K1 was large and K2
% was small. She also assumed K3 was large.

clear
close all

% nK = 50;
% K1 = logspace(-2,4,nK)';
% K2 = logspace(-8,2,nK)';
% K3 = logspace(-2,4,nK)';
% Use the numbers above if you're checking out extremes. Otherwise, the
% ones below.

nK = 50;
K1 = logspace(-2,6,nK)';
K2 = logspace(-6,2,nK)';
K3 = logspace(-4,0,nK)';
F = [2 10 100 1000];
nF = length(F);

%% ========================================================================
% Set up the simulation to show that Lea's simplifications make sense. They
% do.
% =========================================================================
%{
K1 = 100; K2 = 1/K1^2; K3 = K1;
r = 1;
p = [r K1 K2 K3];
F = 1.5;

%
% Get initial conditions
%
ftnhand = @ftn_goentoro_rescale;
tspan = [0 15];
% opts = odeset('reltol',1e-8,'abstol',1e-8);
[t,Y] = ode15s(ftnhand,tspan,[1 1],[],p,1);
figure
plot(t,[Y(:,1) Y(:,2)])
% ylim([0 1.1])

%
% Run the simulation
%
Y0 = Y(end,:)';
[t,Y] = ode15s(ftnhand,tspan,Y0,[],p,F);

figure
plot(t,Y(:,2))
hold on
plot(xlim,Y(1,2)*[1 1],'k :')

%}


%% ========================================================================
% Set up the simulation to show that MY region of param space with NPA is
% boring. It is.
% =========================================================================
%{
K1 = 0.01; K2 = 1; K3 = 1;
r = 1;
p = [r K1 K2 K3];
F = 10;

%
% Get initial conditions
%
ftnhand = @ftn_goentoro;
tspan = [0 15];
% opts = odeset('reltol',1e-8,'abstol',1e-8);
[t,Y] = ode15s(ftnhand,tspan,[1 1],[],p,1);
figure
plot(t,[Y(:,1) Y(:,2)])
% ylim([0 1.1])

%
% Run the simulation
%
Y0 = Y(end,:)';
[t,Y] = ode15s(ftnhand,tspan,Y0,[],p,F);

figure
plot(t,Y(:,2))
hold on
plot(xlim,Y(1,2)*[1 1],'k :')

f = (Y(end,2)-Y(1,2))/Y(1,2);
disp(f)

%}


%% ========================================================================
% Now we look at the behavior of the full model. since we're only
% interested in the steady state values, we don't even need to run a
% simulation. We just calculate the steady states analytically. The
% parameter "r" doesn't matter. We will explore the parameter space to find
% perfect adaptation.
% =========================================================================
% {


%
% setting up for loops to do a parameter sweep
%
K11 = repmat(K1,1,nK); % nK-by-nK
y0 = 1./(K11 + 1); % nK-by-nK
z = zeros(nK,nK,nK,nF);
z0 = zeros(nK,nK,nK);
for i = 1:nK
	
	%
	% Initial condition: F = 1;
	%
	K12 = K1*(1./K2'); 
	z0(:,:,i) = 1./(K11 + 1 + y0.*(K12 + K11/K3(i))); % nK-by-nK
	
	for j = 1:nF
		y = F(j)./(K11 + F(j)); % nK-by-nK
		z(:,:,i,j) = F(j)./(K11 + F(j) + y.*(K12 + F(j)*K11/K3(i)));
		
		
	end
end

%
% f = (z - z0)/z0. And after we calculate this, how to visualize the data?
% Make a threshold.
%
f = (z - repmat(z0,1,1,1,nF))./repmat(z0,1,1,1,nF);
thresh = 0.05;
V = abs(f) < thresh;

%
% Visualize all cases.
% 
for j = 1:nF
	v = V(:,:,:,j);
	
	K11 = repmat(K1,1,nK,nK); K11 = K11(v); K11 = K11(:);
	K22 = repmat(K2',nK,1,nK); K22 = K22(v); K22 = K22(:);
	K33 = permute(K3,[2 3 1]); % rotated K3 into the 3rd dimension
	K33 = repmat(K33,nK,nK,1); K33 = K33(v); K33 = K33(:);
	
	figure
	plot3(K11,K22,K33,'.')
	set(gca,'xscale','log','yscale','log','zscale','log')
	xlabel('K1')
	ylabel('K2')
	zlabel('K3')
	title(['F = ',num2str(F(j))])
	% surf(K1,K2,f(:,:,1,1))
	% set(gca,'xscale','log','yscale','log')
end


%
% Apparently, if K1 is low enough and K2 is high enough, it doesn't matter
% what K3 is. Checking it out here.
% 
figure
S = {'.','o','s','p'};
for j = 1:nF
	v = V(:,:,:,j);
	v = all(v,3);
	
	K11 = repmat(K1,1,nK); K11 = K11(v); K11 = K11(:);
	K22 = repmat(K2',nK,1); K22 = K22(v); K22 = K22(:);
	
	loglog(K11,K22,S{j})
	hold on
	% surf(K1,K2,f(:,:,1,1))
	% set(gca,'xscale','log','yscale','log')
end
xlabel('K1')
ylabel('K2')
title('K3 doesn''t matter')
legend({['F = ',num2str(F(1))],num2str(F(2)),num2str(F(3)),num2str(F(4))},...
	'location','southeast')

%}


%% ========================================================================
% Since there is a very finely-tuned surface in param space where near-
% perfect adaptation occurs, we will use design rules to see if that is the
% surface we're finding.
%
% UPDATE: Turns out, yes it is.
% =========================================================================
%{

%
% Below is the design rule. We can find the K3 that solves this equation
% for all values of K1 and K2.
%
% 0 = K1 + K1*F*(1/K2*(1/(K1+F) - 1/(K1+1)) + ...
%	1/K3*(F/(K1+F) - 1/(K1+1)) - 1)
%
% We will loop through all values of F
%
K332 = K1*K2';
for i = 1:nF
	F1 = F(i);
	
	%
	% Solving for K3
	%
	LHS = (1-1/F1) - (1./(K1 + F1) - 1./(K1 + 1))*(1./K2');
	RHS = F1./(K1 + F1) - 1./(K1 + 1);
	K33 = repmat(RHS,1,nK)./LHS;
	
	figure
% 	loglog(K2,K33')
	surf(K1,K2,K33')
% 	hold on
% 	surf(K1,K2,K332')
	set(gca,'xscale','log','yscale','log','zscale','log')
	xlabel('K_1')
	ylabel('K_2')
	zlabel('K_3*')
	title(['F = ',num2str(F1)])

% 	%
% 	% Solving for K2 instead. Why doesn't this work as well? Oh, it's
% 	% because the surface for K2 goes almost straight up, and for K3 large
% 	% enough, there is no solution. In fact, if you look at the "surf" plot
% 	% of K33 vs K1,K2, K33 never gets above 1.0, so there will be a big
% 	% region of the K1,K3 plane where K22 won't have a solution. It's a
% 	% vertical asymptote? Then, since this is on a log scale, the negative
% 	% data on the other side of that asymptote are not plotted.
% 	%
% 	LHS = (1-1/F1) - (F1./(K1 + F1) - 1./(K1 + 1))*(1./K3');
% 	RHS = 1./(K1 + F1) - 1./(K1 + 1);
% 	K22 = repmat(RHS,1,nK)./LHS;
% 	
% 	figure
% 	surf(K1,K3,K22')
% 	set(gca,'xscale','log','yscale','log','zscale','log')
% 	xlabel('K_1')
% 	ylabel('K_3')
% 	zlabel('K_2*')
% 	title(['F = ',num2str(F1)])
end


%}

%% ========================================================================
% Alternate of the above: if K3 = K1*K2, then we reduce the dimensionality.
% 
% UPDATE: Actually, this doesn't work. If K3 = K1*K2 (ie, if there's no
% cooperativity), then it's impossible to get perfect adaptation. Weird.
% =========================================================================
%{

%
% K2*(1 - 1/F) = (1/K1)*(F/(K1+F) - 1/(K1+1)) + (1/(K1+F) - 1/(K1+1))
%
% We will loop through all values of F
%
figure
for i = 1:nF
	F1 = F(i);
	LHS = (1-1/F1);
	RHS = (1./K1).*(F1./(K1 + F1) - 1./(K1 + 1)) + (1./(K1 + F1) - 1./(K1 + 1));
	K22 = RHS./LHS;
	
	semilogx(K1,K22)
	hold on
	xlabel('K_1')
	ylabel('K_2')
end
legend({['F = ',num2str(F(1))],num2str(F(2)),num2str(F(3)),num2str(F(4))},...
	'location','southeast')



%}

%% ========================================================================
% We have confirmed that the very finely-tuned surface in param space where
% near-perfect adaptation occurs is indeed coincident with the design
% rules. Now we will enact those design rules, tweak them a bit, and see if
% we still have good adaptation.
%
% In this cell, we are varying K1 and K2 across many logs in the parameter
% space (see cell 1), and forcing K3 to the design rule for perfect
% adaptation. Then, we double one of the parameters and calculate f, which
% is defined as (z_final - z_initial)/z_initial. Note that our starting
% point, before we double one of the parameters, has been tuned so that f
% equals zero. So any f that we calculate here will also be a delta_f.
%
% UPDATE: There is a weird region of parameter space that is robust to all
% paramter doublings. To see this, look at the "max intensity" projection
% of all three paramter doublings (K1,K2,K3) and the F1 being 10*F1. Only
% low K1 (0.05 -- 0.15, depending on K2) and a corner of the param space of
% very high K1 (greater than 100) and super-low K2 (less than 1e-4) survive
% all paramter doublings.
%
% ANOTHER UPDATE: this cell is slightly obsoleted by the next cell, which
% does the same thing, with more updated standardized code, but the next
% cell also looks at the neg fbk loop.
% =========================================================================
%{

F0 = 1;
F1 = 10;
LHS = (1-1/F1) - (1./(K1 + F1) - 1./(K1 + 1))*(1./K2');
RHS = F1./(K1 + F1) - 1./(K1 + 1);
K331 = repmat(RHS,1,nK)./LHS;


%
% Now we run the for loop. The first time through, it will be baseline. The
% second through third times, we will double one of the constants.
%
exprn = {'''all normal'';','K111 = 2*K1;','K222 = 2*K2;','K33 = 2*K331;','F11 = 2*F1;'};
f_out = zeros(nK,nK);
for i = 1:length(exprn)
	K111 = K1; K222 = K2; K33 = K331; F11 = F1;
	eval(exprn{i})

	%
	% setting up for loops to do a parameter sweep
	%
	K11 = repmat(K111,1,nK); % nK-by-nK
	K22 = repmat(K222,1,nK);
	
	%
	% Initial condition: F = 1;
	%
	y0 = F0./(K11 + F0); % nK-by-nK
	K12 = K111*(1./K222'); % nK-by-nK
	z0 = 1./(K11 + 1 + y0.*(K12 + K11./K33)); % nK-by-nK
	
	y = F11./(K11 + F11); % nK-by-nK
	z = F11./(K11 + F11 + y.*(K12 + F11*K11./K33));
	
	
	%
	% f = (z - z0)/z0. And after we calculate this, how to visualize the data?
	% Make a threshold.
	%
	f = (z - z0)./z0;
	thresh = 0.05;
	thresh2 = 0.2;
	V = abs(f) < thresh;
	f_out = max(f_out,f);
	
% 	%
% 	% plotting
% 	%
% 	figure
% 	f(1,1) = 0; if i > 1, f(end,end) = 0.75; end
% 	pcolor(K1,K2,abs(f'))
% 	set(gca,'xscale','log','yscale','log')
% 	if i > 1
% 		shading interp
% 	else
% 		shading flat
% 	end
% 	colorbar
% 	xlabel('K_1')
% 	ylabel('K_2')
% 	title(exprn{i})
% 	hold on
% 	
% 	c = contourc(K1,K2,abs(f'),[thresh thresh2]);
% 	if ~isempty(c)
% 		[X,Y,C] = extractcontour(c);
% 		for j = 1:length(Y)
% 			if C(j) == thresh
% 				plot(X{j},Y{j},'w')
% 			elseif C(j) == thresh2
% 				plot(X{j},Y{j},'m')
% 			end
% 		end
% 	end

end

%
% plotting "max intensity projection"
%
figure
f_out(1,1) = 0; f_out(end,end) = 0.75;
pcolor(K1,K2,abs(f_out'))
set(gca,'xscale','log','yscale','log')
shading interp
colorbar
xlabel('K_1')
ylabel('K_2')
title('max')
hold on

c = contourc(K1,K2,abs(f_out'),[thresh thresh2]);
if ~isempty(c)
	[X,Y,C] = extractcontour(c);
	for j = 1:length(Y)
		if C(j) == thresh
			plot(X{j},Y{j},'w')
		elseif C(j) == thresh2
			plot(X{j},Y{j},'m')
		end
	end
end



%}



%% ========================================================================
% So, not only do we have a finely-tuned surface in the K[1-3] parameter
% space, but we also have a limited fraction of that space in which the
% adaptation property is maintained upon perturbation to the parameters. So
% can the addition of negative feedback generally rescue this fragility?
% Here we will test that question.
%
% To do this, we will sweep the parameter space in K1,K2,K3, while allowing
% K4 and K8 to take on low, med, and high values. We will fix F = 10. Since
% we are only looking at the steady states (initially, here), we don't care
% about r or s, so they will generically be set to 1.
% =========================================================================
%{

nw = 1; % number of variations in K4 and K8.
% K4 = logspace(-2,0,nw)'; 
K4 = 1;
K8 = K4;

%
% setting up for loops to do a parameter sweep
%
F0 = 1;
F1 = 10;
z = zeros(nK,nK,nK,nw^2);
z0 = zeros(nK,nK,nK,nw^2);

fhandle = @ftn_goentoro_wfbk_SS;
zguess = 0.5;
for i = 1:nK
	K11 = K1(i);
	y0 = F0/(K11 + F0);
	y1 = F1/(K11 + F1);
	A0 = F0/K11;
	A1 = F1/K11;
	
	for j = 1:nK
		K22 = K2(j);
		for k = 1:nK
			K33 = K3(k);
			
			B0 = 1 + F0/K11 + y0/K22 + F0.*y0/K33;
			B1 = 1 + F1/K11 + y1/K22 + F1.*y1/K33;
			
			count = 1;
			for o = 1:nw
				K44 = K4(o);
				K5 = K11*K44; K6 = K22*K44; K7 = K33*K44;
				
				C0 = 1/K44 + F0/K5 + y0/K6 + F0.*y0/K7;
				C1 = 1/K44 + F1/K5 + y1/K6 + F1.*y1/K7;
			
				for p = 1:nw
					K88 = K8(p);
					
					
					z0(i,j,k,count) = -(K88.*B0-A0)/2./(B0+C0) + ...
						sqrt((K88.*B0-A0).^2 + 4*(B0+C0).*A0.*K88)/2./(B0+C0);
					z(i,j,k,count) = -(K88.*B1-A1)/2./(B1+C1) + ...
						sqrt((K88.*B1-A1).^2 + 4*(B1+C1).*A1.*K88)/2./(B1+C1);
					
					count = count + 1;
				end				
			end	
		end
	end
	disp(i)
end

%
% f = (z - z0)/z0. And after we calculate this, how to visualize the data?
% Make a threshold.
%
f = (z - z0)./z0;
thresh = 0.05;
V = abs(f) < thresh;

%
% Visualize all cases.
% 
figure
count = 1;
for i = 1:nw
	for j = 1:nw
	v = V(:,:,:,count);
	
	K11 = repmat(K1,1,nK,nK); K11 = K11(v); K11 = K11(:);
	K22 = repmat(K2',nK,1,nK); K22 = K22(v); K22 = K22(:);
	K33 = permute(K3,[2 3 1]); % rotated K3 into the 3rd dimension
	K33 = repmat(K33,nK,nK,1); K33 = K33(v); K33 = K33(:);
	
	subplot(nw,nw,count)
	plot3(K11,K22,K33,'.')
	set(gca,'xscale','log','yscale','log','zscale','log')
	xlabel('K1')
	ylabel('K2')
	zlabel('K3')
	title(['K4 = ',num2str(K4(i)),'; K8 = ',num2str(K8(j))])
	% surf(K1,K2,f(:,:,1,1))
	% set(gca,'xscale','log','yscale','log')
	
		count = count + 1;
	end
end



%}



%% ========================================================================
% After doing some more hand-written calculations, it turns out the
% presence of the neg fbk loop does not change the PA-surface at all. IOW,
% if I vary K1 and K2, I can calculate a surface of K3 that achieves PA.
% (ie, a design constraint). This surface doesn't change at all with the
% addition of the NFL.
%
% The purpose of this cell is to co-plot with the previous cell. So have
% them both on if you want that.
% =========================================================================
%{

F0 = 1;
F1 = 10;
LHS = (1/F0-1/F1) - (1./(K1 + F1) - 1./(K1 + F0))*(1./K2');
RHS = F1./(K1 + F1) - F0./(K1 + F0);
K33 = repmat(RHS,1,nK)./LHS;


hold on
surf(K1,K2,K33')
set(gca,'xscale','log','yscale','log','zscale','log')
xlabel('K_1')
ylabel('K_2')
zlabel('K_3*')
title(['F = ',num2str(F1)])

%}





%% ========================================================================
% Now that we know that the K3* surface is the same, we can look at the
% gradient in "f" to see if NPA is more robust w/NFL vs w/o
%
% NOTE: this is the old version where we calculate the derivative
% analytically, but now we are lazy and calculate it by fin diff. We can do
% this for all three K's. Also, we added back an epsln.
% =========================================================================
%{

epsln = 0.05;
F0 = 1;
F1 = 10;
LHS = (1/F0-1/F1) - (1./(K1 + F1) - 1./(K1 + F0))*(1./K2');
RHS = F1./(K1 + F1) - F0./(K1 + F0);
K33 = repmat(RHS,1,nK)./LHS;

K11 = repmat(K1',nK,1);
K22 = repmat(K2,1,nK);

y0 = F0./(K11 + F0);
y1 = F1./(K11 + F1);


nw = 4; % number of variations in K4 and K8.
K4 = logspace(-3,0,nw)'; 
K8 = K4;

% K4 = 1;
% K8 = 1;
% K5 = K11*K4; K6 = K22*K4; K7 = K33*K4;


%
% Without the neg fbk loop
%
%
z0 = F0./K11./(1 + F0./K11 + y0./K22 + F0.*y0./K33); % nK-by-nK
dz0dK3 = K11.*y0.*z0.^2./K33.^2;

z1 = F1./K11./(1 + F1./K11 + y1./K22 + F1.*y1./K33);
dz1dK3 = K11.*y1.*z1.^2./K33.^2;

dfOLdK3 = (dz1dK3.*z0 - z1.*dz0dK3)./z0^2;
DK_OL = epsln./abs(dfOLdK3);
rho_OL = 1./(K33.*abs(dfOLdK3));


%
% For loops to change K4,K8 to determine if neg fbk helps, depending on
% strength of neg fbk.
%
Dmax = 0;
Dmin = 0;
D_out = cell(nw,nw);
DK_CL = cell(nw,nw);
rho_out = cell(nw,nw);
rho_CL = cell(nw,nw);
for i = 1:nw
	K44 = K4(i);
	K5 = K11*K44; K6 = K22*K44; K7 = K33*K44;
	
	for j = 1:nw
		
		%
		% With the neg fbk loop
		%
		% initial state
		K88 = K8(j);
		A0 = F0./K11;
		B0 = 1 + A0 + y0./K22 + F0.*y0./K33;
		C0 = 1/K44 + F0./K5 + y0./K6 + F0.*y0./K44./K33;
		z0 = -(K88.*B0-A0)/2./(B0+C0) + ...
			sqrt((K88.*B0-A0).^2 + 4*(B0+C0).*A0.*K88)/2./(B0+C0);
		dB0dK3 = -F0.*y0./K33.^2;
		dC0dK3 = -F0.*y0./K33.^2/K44;
		dz0dK3 = ((dB0dK3 + dC0dK3).*z0.^2 + K88*dB0dK3.*z0) ./ ...
			(-(K88.*B0 - A0) - 2*(B0 + C0).*z0);
		
		% final state
		A1 = F1./K11;
		B1 = 1 + A1 + y1./K22 + F1.*y1./K33;
		C1 = 1/K44 + F1./K5 + y1./K6 + F1.*y1./K44./K33;
		z1 = -(K88.*B1-A1)/2./(B1+C1) + ...
			sqrt((K88.*B1-A1).^2 + 4*(B1+C1).*A1.*K88)/2./(B1+C1);
		dB1dK3 = -F1.*y1./K33.^2;
		dC1dK3 = -F1.*y1./K33.^2/K44;
		dz1dK3 = ((dB1dK3 + dC1dK3).*z1.^2 + K88*dB1dK3.*z1) ./ ...
			(-(K88.*B1 - A1) - 2*(B1 + C1).*z1);
		
		dfCLdK3 = (dz1dK3.*z0 - z1.*dz0dK3)./z0^2;
		
		DK_CL{i,j} = epsln./abs(dfCLdK3);
		rho_CL{i,j} = 1./(K33.*abs(dfCLdK3));
		
		%
		% Relative difference in allowable change in K3
		%
% 		D = abs(dfOLdK3) - abs(dfCLdK3);
% 		Dmax = max(max(D(:)),Dmax);
% 		Dmin = min(min(D(:)),Dmin);
% 		D_out{i,j} = D;

		
		D = (DK_CL{i,j} - DK_OL)./K33;
		Dmax = max(max(D(:)),Dmax);
		Dmin = min(min(D(:)),Dmin);
		D_out{i,j} = D;
		
		
		rho_out{i,j} = rho_CL{i,j} - rho_OL;
		
	end
end

Dtop = 100;
COLRS = {'y' 'w' 'r' 'b'};
thresh = [0 5 10 100];

%
% plotting the OL robustness
%
% D = DK_OL./K33;
D = rho_OL;
c = contourc(K1,K2,D,thresh);
D(D > Dtop) = Dtop;

figure
D(end,end) = 0;
pcolor(K1,K2,D)
set(gca,'xscale','log','yscale','log')
shading interp
colorbar
xlabel('K_1')
ylabel('K_2')
title('open loop robustness')
hold on

% contours
if ~isempty(c)
	[X,Y,C] = extractcontour(c);
	for k = 1:length(Y)
		vthresh = C(k) == thresh;
		plot(X{k},Y{k},COLRS{vthresh})		
			
% 		if C(k) == thresh
% 			plot(X{k},Y{k},'w')
% 		elseif C(k) == thresh2
% 			plot(X{k},Y{k},'m')
% 		end
	end
end



%
% plotting the CL robustness
%
figure
count = 1;
for i = 1:nw
	
	for j = 1:nw
		
% 		D = DK_CL{i,j}./K33;
		D = rho_CL{i,j};
		c = contourc(K1,K2,D,thresh);
		D(end,end) = 0;
		D(D < 0) = 0;
		D(D > Dtop) = Dtop;
		
		subplot(nw,nw,count);
% 		D(1,1) = Dmin;
% 		D(end,1) = Dmax;
		pcolor(K1,K2,D)
		set(gca,'xscale','log','yscale','log')
		shading interp
		colorbar
		xlabel('K_1')
		ylabel('K_2')
		title(['K4 = ',num2str(K4(i)),'; K8 = ',num2str(K8(j))])
		hold on
		
		
		if ~isempty(c)
			[X,Y,C] = extractcontour(c);
			for k = 1:length(Y)
				vthresh = C(k) == thresh;
				plot(X{k},Y{k},COLRS{vthresh})
			end
		end
		
		count = count + 1;
	end
end


%
% plotting the difference
%
figure
count = 1;	
for i = 1:nw
	
	for j = 1:nw
		
		D = rho_out{i,j};
		c = contourc(K1,K2,D,thresh);
		
		D(end,end) = 0;
		D(D < 0) = 0;
		D(D > Dtop) = Dtop;
		
		subplot(nw,nw,count);
% 		D(1,1) = Dmin;
% 		D(end,1) = Dmax;
		pcolor(K1,K2,D)
		set(gca,'xscale','log','yscale','log')
		shading interp
		colorbar
		xlabel('K_1')
		ylabel('K_2')
		title(['K4 = ',num2str(K4(i)),'; K8 = ',num2str(K8(j))])
		hold on
		
		
		if ~isempty(c)
			[X,Y,C] = extractcontour(c);
			for k = 1:length(Y)
				vthresh = C(k) == thresh;
				plot(X{k},Y{k},COLRS{vthresh})
			end
		end
		
		count = count + 1;
	end
end

%}