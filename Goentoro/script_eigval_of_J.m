% script_eigval_of_J
%
% This script attempts to see how the value of J changes the eigenvalues of
% our system.

clear
close all

% F1 = 10;
% K1 = 1; 
% K2 = 0.3;
% K4 = 0.1;
% K8 = 0.1;
% n = 2;

%
% In this first section, we acquire the beta that gives us the instability
% point. This "betastar" is the frequency of oscillations at instability.
% Note that this beta (and subsequently, J) is ONLY dependent on tauw,
% tauz, thetaw, and thetaz. That is an interesting observation, because the
% steady states are independent of those parameters.
%
tauw = 1;
tauz = 1;
thetaz = 0.5;
thetaw = 0.5;

betamax = pi/2/(thetaw + thetaz) - 0.05;
beta = linspace(0,betamax)';
fbeta = eigval_alpha0(beta,tauz,tauw,thetaz,thetaw);

betastar = advNtnDU(@eigval_alpha0,betamax/2,[1e-4,betamax],[],[],[],[],[],...
	tauz,tauw,thetaz,thetaw);

%
% Now that we have beta, we calculate J, which is the other varible
% parameter in the characteristic equation. In a very complicated way, all
% other parameters (K1...K8, as well as the steady state values of our
% variables) are embedded in J.
%
Jstar = (1-tauw*tauz*betastar.^2)./cos((thetaw+thetaz).*betastar);

%% ========================================================================
% Running a for loop (w/dumb continuation) to trace out the eigenvalues as
% a function of J. We do this twice: once in one direction, another time in
% the other direction.
% =========================================================================

nJ = 250;
J = [linspace(Jstar,0,nJ)' linspace(Jstar,-100,nJ)'];
Alpha = zeros(nJ,2);
Beta = zeros(nJ,2);
tolerance = 1e-6; nStepsMax = 100;
delt = 1e-4;
for k = 1:2
	dJ = J(2,k) - J(1,k);
	X = [0; betastar];
		
	%
	% Now we run the for loop. The Newton while loop is embedded.
	%
	for i = 1:nJ
		
		dX = 1; F = 0; nSteps = 1;
		while norm([dX;F]) > tolerance && nSteps < nStepsMax
			% 	X = fsolve(@ftn_chareqn,X0,[],tauz,tauw,thetaz,thetaw,J(i));
			alpha = X(1);
			beta = X(2);
			
			%
			% Calculate F
			%
			F = ftn_chareqn(X,tauz,tauw,thetaz,thetaw,J(i,k));
			
			%
			% Calc "Jac"
			%
			alpha1 = alpha + delt;
			F1 = ftn_chareqn([alpha1;beta],tauz,tauw,thetaz,thetaw,J(i,k));
			F_alpha = (F1 - F)/delt;
			
			beta1 = beta + delt;
			F1 = ftn_chareqn([alpha;beta1],tauz,tauw,thetaz,thetaw,J(i,k));
			F_beta = (F1 - F)/delt;
			Jac = [F_alpha F_beta];
			
			%
			% Take Newton step
			%
			dX = -Jac \ F;
			X = X + dX;
			nSteps = nSteps + 1;
			
		end
		
		Alpha(i,k) = X(1);
		Beta(i,k) = X(2);
		
		%
		% Calc "G" at Jstar
		%
		J1 = J(i,k) + delt;
		F1 = ftn_chareqn(X,tauz,tauw,thetaz,thetaw,J1);
		G = (F1 - F)/delt;
		
		%
		% Predictor
		%
		dXdJ = -Jac \ G;
		X = X + dXdJ*dJ;
		
	end

end

%
% Now that we have these, let's combine them together
%
J = J(:);
Alpha(1) = NaN;
Alpha = Alpha(:);
Beta = Beta(:);
[J,isort] = sort(J);
Alpha = Alpha(isort);
Beta = Beta(isort);
v = isnan(Alpha);
Alpha(v) = [];
Beta(v) = [];
J(v) = [];

figure
plot(Alpha,Beta)

save Mat/2018-11-30_eigval_of_J J Alpha Beta betastar Jstar tauw tauz thetaw thetaz










