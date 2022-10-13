% script_regionsI_and_II
%
% This script finds the boundaries of Regions I and II in the K1,K2 plane.

delt = 1e-4;

%% ========================================================================
% First we must attempt to delineate the boundaries of Regions I and II.
% For Region I, the boundary is evaluated at K12 -> +infty, but you have to
% solve implicitly for K1 vs K2, using Newton's method. However, the method
% is having a really hard time solving for all points, so we'll start at an
% "arbitrary" place, which we happen to know works, and it's right on the
% boundary for where we want to be, and we'll do poor-man's continuation
% forward in K2, but we have to use a more dense mesh than before.
% =========================================================================

% {

nK_1 = 300;
K2_I_FBL = logspace(log10(0.1485),2,nK_1)';
K1_I_FBL = NaN(nK_1,1);
p = [1 1 1 Inf K3 K4 1 1];
x = log(0.0102);
for i = 1:nK_1
	K2 = K2_I_FBL(i);
	p(3) = K2;
	
	y0 = @(K1)x0./(K1 + x0);
	y1 = @(K1)x1./(K1 + x1);
	fhandle = @(K1,K2,A)(A/(1 + vareps)*(K1.*y0(K1)/x0) - (K1.*y1(K1)/x1)) / ...
		((1 + K1/x1) - A/(1 + vareps)*(1 + K1/x0)) - K2;
	
	f = 1; dx = 0; nSteps = 0;
	tolerance = 1e-6; nStepsMax = 500;
	while norm([f;dx]) > tolerance && nSteps < nStepsMax
		%
		% calc f
		%
		p(2) = exp(x);
		[~,~,w0] = ftn_goentoro_ss(p,x0,n);
		[~,~,w1] = ftn_goentoro_ss(p,x1,n);
		A = (1 + w0/K3)./(1 + w1/K3);
		f = fhandle(exp(x),K2,A);
		
		%
		% Calc fprime
		%
		x_ = x + delt;
		p(2) = exp(x_);
		[~,~,w0] = ftn_goentoro_ss(p,x0,n);
		[~,~,w1] = ftn_goentoro_ss(p,x1,n);
		A = (1 + w0/K3)./(1 + w1/K3);
		f1 = fhandle(exp(x_),K2,A);
		fprime = (f1 - f)/delt;
		
		%
		% Take the step
		%
		dx = -f/fprime;
		x = x + dx;
		nSteps = nSteps + 1;
	end
	
	if nSteps < nStepsMax
		K1_I_FBL(i) = exp(x);
	end
	
	if i < nK_1%142 % serendipitous that this only works for index 142 or lower
		
		%
		% Now calc. the derivative of f wrt K2
		%
		K2_1 = K2*(1 + delt);
		p(3) = K2_1;
		[~,~,w0] = ftn_goentoro_ss(p,x0,n);
		[~,~,w1] = ftn_goentoro_ss(p,x1,n);
		A = (1 + w0/K3)./(1 + w1/K3);
		f1 = fhandle(exp(x),K2_1,A);
		G = (f1 - f)/delt;
		dlogK1dlogK2 = - G/fprime;
		dlogK2 = log(K2_I_FBL(i+1)) - log(K2);
		x = x + dlogK1dlogK2*dlogK2; % next guess
	end
	1;
end

%}

%% ========================================================================
% Now Region II. Here, the boundary is also evaluated at K12 -> +infty, and
% again you have to solve implicitly for K1 vs K2, using Newton's method.
% However, the method is having a really hard time solving for all points,
% so we'll start at an "arbitrary" place, which we happen to know works,
% and it's right on the boundary for where we want to be, and we'll do
% palc, starting backward in K1.
% =========================================================================
% {
nK_1 = 300;
K1_II_FBL = zeros(nK_1,1);
K2_II_FBL = zeros(nK_1,1);
p = [1 1 1 Inf K3 K4 1 1];
K1 = 5.5943e+04;
K2 = 1.005071423775046e-06;
x = log([K1;K2]);
ds = 0.1;

y0 = @(K1)x0./(K1 + x0);
y1 = @(K1)x1./(K1 + x1);
fhandle = @(K1,K2,A)A/(1 + vareps)*(K2 + K1*K2/x0 + K1.*y0(K1)/x0) - ...
	(K2 + K1*K2/x1 + K1.*y1(K1)/x1);
for i = 1:nK_1
	
	F = 1; dx = 0; nSteps = 0;
	tolerance = 1e-6; nStepsMax = 50;
	while norm([F;dx]) > tolerance && nSteps < nStepsMax
		
		%
		% calc f
		%
		K1 = exp(x(1)); K2 = exp(x(2));
		p(2) = K1; p(3) = K2;
		[~,~,w0] = ftn_goentoro_ss(p,x0,n);
		[~,~,w1] = ftn_goentoro_ss(p,x1,n);
		A = (1 + w0/K3)./(1 + w1/K3);
		f = fhandle(K1,K2,A);
		if i > 1
			N = dxds0'*(x - x_0) - ds;
			F = [f;N];
		else
			F = f;
		end
		
		%
		% Now calc. the derivative of f wrt x
		%
		dfdx = [0 0];
		for k = 1:2
			x_ = x;
			x_(k) = x_(k) + delt;
			p1 = p;
			p1(k+1) = exp(x_(k));
			[~,~,w0] = ftn_goentoro_ss(p1,x0,n);
			[~,~,w1] = ftn_goentoro_ss(p1,x1,n);
			A = (1 + w0/K3)./(1 + w1/K3);
			f1 = fhandle(exp(x_(1)),exp(x_(2)),A);
			dfdx(k) = (f1 - f)/delt;
		end
		% 		if i == 1
		% 			dfdx(2) = [];
		% 		end
		
		%
		% Derivatives of N wrt K1 and K2
		%
		if i > 1
			dNdx = dxds0;
		end
		
		%
		% Put together the Jacobian
		%
		if i == 1
			J = dfdx(1);
		else
			J = [dfdx;dNdx'];
		end
		
		%
		% Take the step
		%
		dx = -J \ F;
		if i == 1, dx1 = [dx;0]; dx = dx1; end
		x = x + dx;
		nSteps = nSteps + 1;
	end
	
	if nSteps < nStepsMax
		K1_II_FBL(i) = exp(x(1));
		K2_II_FBL(i) = exp(x(2));
	end
	
	if i < nK_1
		
		x_0 = x;
		dx1dx2 = -dfdx(1) \ dfdx(2);
		
		if i == 1
			%
			% Initial calculation of dx2ds
			%
			dx2ds = 1/(sqrt(1 + dx1dx2^2));
		else
			%
			% Typical...watch out for sign
			%
			dx2ds = 1/(dxds0(1)*dx1dx2 + dxds0(2));
		end
		dx1ds = dx1dx2*dx2ds;
		dxds0 = [dx1ds; dx2ds];
		
		
		x = x_0 + dxds0*ds;
	end
	1;
end

save Mat/regionsI_and_II

