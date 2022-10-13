% script_regionsI_II_coop
%
% This script finds the boundaries of Regions I and II in the K1,K2 plane.
% This is for different cooperativity levels of X with W and Y with W.


%% ========================================================================
% First we must attempt to delineate the boundaries of Regions I and II.
% For Region I, the boundary is evaluated at K12 -> +infty, but you have to
% solve implicitly for K1 vs K2, using Newton's method. However, the method
% is having a really hard time solving for all points, because there is a
% K1_upper that is messing us up. So we'll have to solve for that K1_upper
% first, and use that as an upper bound for advNtn. Then we'll do
% poor-man's continuation forward in K2, but we have to use a more dense
% mesh than before.
% =========================================================================

% {

%
% Solve for K1_upper. Region I cannot have a K1 greater than this. Another
% way to think about it is that Region I has a vertical asymptote at this
% value of K1.
%
fhandle = @ftn_RegionI_K1upper_w_coop;

% % 
% % See what the landscape of "f" looks like
% % 
% x_try = linspace(-10,10,150);
% f_try = zeros(length(x_try),1);
% for i_try = 1:length(x_try)
% 	f_try(i_try) = fhandle(x_try(i_try),K3,K4,C13,C23,x0,x1,n,vareps);
% end
% figure
% semilogx(exp(x_try),f_try)

x = 0;
bounds = [-10 10];
nStepsMax = 100;
[xupper,fprime,f,nSteps] = ...
	advNtnDU(fhandle,x,bounds,[],nStepsMax,[],[],[],...
	K3,K4,C13,C23,x0,x1,n,vareps);
K1_upper = exp(xupper);

%
% For loop for continuation
%
nK_1 = 300;
K2_I_FBL = logspace(-6,2,nK_1)';
K2_I_FBL = flipud(K2_I_FBL);
x = xupper*(1 - sign(xupper)*delt);
K1_I_FBL = NaN(nK_1,1);
fhandle = @ftn_RegionI_II_FB_w_coop;
for i = 1:nK_1
	K2 = K2_I_FBL(i);
	
% 	%
% 	% See what the landscape of "f" looks like
% 	%
% 	x_try = linspace(-5,5,150);%xupper*(1+delt),150);
% 	f_try = zeros(length(x_try),1);
% 	for i_try = 1:length(x_try)
% 		f_try(i_try) = fhandle(x_try(i_try),K2,K3,K4,C13,C23,x0,x1,n,vareps);
% 	end
% 	figure
% 	semilogx(exp(x_try),f_try)
	
	
	%
	% Do advntn
	%
	bounds = [-10 xupper];
	[x,fprime,f,nSteps] = ...
		advNtnDU(fhandle,x,bounds,[],nStepsMax,[],[],[],...
		K2,K3,K4,C13,C23,x0,x1,n,vareps);
	if nSteps < nStepsMax
		K1_I_FBL(i) = exp(x);
	end
	if K1_I_FBL(i) < 1e-2
		break
	end
	
	if i < nK_1%142 % serendipitous that this only works for index 142 or lower
		
		%
		% Now calc. the derivative of f wrt K2
		%
		K2_1 = K2*(1 + delt);
		f1 = fhandle(x,K2_1,K3,K4,C13,C23,x0,x1,n,vareps);
		G = (f1 - f)/delt;
		dlogK1dlogK2 = - G/fprime;
		dlogK2 = log(K2_I_FBL(i+1)) - log(K2);
		x = x + dlogK1dlogK2*dlogK2; % next guess
	end
	1;
end
K2_I_FBL(i+1:end) = [];
K1_I_FBL(i+1:end) = [];
K2_I_FBL = flipud(K2_I_FBL);
K1_I_FBL = flipud(K1_I_FBL);



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
p = [1 1 1 Inf K3 K4 1 C13 C23];
K1 = 5.5943e+04;
K2 = 1.005071423775046e-06;
x = log([K1;K2]);
ds = 0.1;

for i = 1:nK_1
	
	F = 1; dx = 0; nSteps = 0;
	tolerance = 1e-6; nStepsMax = 50;
	while norm([F;dx]) > tolerance && nSteps < nStepsMax
		
		%
		% calc f
		%
		K1 = exp(x(1)); K2 = exp(x(2));
		p(2) = K1; p(3) = K2;
		
		f = fhandle(x(1),K2,K3,K4,C13,C23,x0,x1,n,vareps);
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
			x_1 = x;
			x_1(k) = x_1(k) + delt;
			p1 = p;
			p1(k+1) = exp(x_1(k));
			
			if k == 1
				f1 = fhandle(x_1(1),K2,K3,K4,C13,C23,x0,x1,n,vareps);
			else
				f1 = fhandle(x(1),p1(3),K3,K4,C13,C23,x0,x1,n,vareps);
			end
			dfdx(k) = (f1 - f)/delt;
		end
		
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
	
	%
	% The PALC predictor step
	%
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

%}