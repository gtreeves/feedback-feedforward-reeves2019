% script_Fig3_FBL_K12NPA_K1K2sweep
%
% This script goes through each of the values of K1,K2 and finds the value
% of K12 that gives NPA. We use a double-for loop and a Newton in the
% middle. We will choose moderate levels of neg fbk strength.

delt = 1e-4;

%% ========================================================================
% Make a for loop to solve, by Ntn's method, the value of K12 that gives f =
% vareps in the neg fbk loop scenario, for many different K1,K2's. Given
% the strange-looking K1,K2 plane, with Regions I and II carved out, we
% will have to construct this plane by hand. We will do this in 2 phases:
% (1a) Start in the lower left corner and do continuation up K1 until we
%	hit Region II. Then move up by one increment of K2 and do it again. We
%	do this until moving up by one increment of K2 (at lowest K1) causes us
%	to hit Region I. NOTE: At some point, we will no longer hit Region II,
%	and instead will continue all the way until we hit max K1.
% (1b) Start back up again with that same value of K2 that we left off
%	with. Increment up by one value of K2. Then do continuation down K1
%	until we hit the min value of K1. Then move up by one increment of K2
%	and do it again. We do this until we hit the max value of K2. NOTE: At
%	some point, we will hit Region I instead of being able to continue all
%	the way until we hit min K1.
% (2) Start back up again with that same value of K2 that we left off with
%	from the end of Phase (1). Increment down by one value of K2. Then do
%	continuation down K1 until we hit Region II. Then move down by one
%	increment of K2 and do it again. We do this until we hit the min value
%	of K2.
% =========================================================================
K12_NPA_FBL = NaN(nK,nK,2);

%
% Need max of and index of max of K2 along Region II.
%
[K2_II_FBL_max,iII_max] = max(K2_II_FBL);
iII = find(K22 < K2_II_FBL_max);
iII = iII(end);

%
% Need index of min of K2 along Region I.
%
iI = find(K22 < K2_I_FBL(1));
iI = iI(end);

%
% OK, this set of loops is a bit crazy and getting out of hand. There are
% two loops over k: one for vareps > 0, one for vareps < 0. Within that, we
% loop over both K1 (inner) and K2 (outer). However, the loop on K2 has two
% Phases: one going up, the other coming down. To make matters more
% complicated, in the middle of Phase 1, the direction that we continue in
% K1 flips (starts out low to high, then switches to high to low). In
% addition, because of the no-go regions, there are several checks to make
% sure we exit the loops (or go to next i, or something) at the right time.
% And those checks differ depending on the Phase and the direction of K1
% movement and whether vareps > 0 or not. When I started writing this loop,
% I thought it would be more "elegant" (ie, readable) if everything was
% contained in one block, with if statements controlling the checks and
% switches. Now I am not so sure...
%
vareps1 = vareps;
for k = 1:2
	if k == 2
		vareps1 = -vareps1;
		
		%
		% If vareps < 0, then there is a lower limit to K1, below which K12 can go
		% all the way down to zero and you still get NPA. IOW, it doesn't make
		% sense to find the lower bound on K12 for NPA if K1 is too low. So we'll
		% skip that part.
		%
		K1_lower = vareps1/(1/x1 - 1/x0*(1 + vareps1));
		
		%
		% Establish this K1_lower as a lower bound.
		%
		jIII = find(K11 <= K1_lower);
		jIII = jIII(end) + 1;
		
		K12_NPA_FBL(:,1:jIII-1,k) = 0;
		
	else
		jIII = 1;		
	end
	
	%
	% Define function handles.
	%
	y0 = @(K1)x0./(K1 + x0);
	y1 = @(K1)x1./(K1 + x1);
	fhandle = @(K1,K2,K12,A)K1*(y1(K1) - A*y0(K1)/(1+vareps1)) / ...
		(A/(1+vareps1)*(1+K1/x0+K1*y0(K1)/K2/x0) - ...
		(1+K1/x1+K1*y1(K1)/K2/x1)) - K12;
	
	%
	% Define jmesh. This goes from 1 to nK, unless you have vareps < 0 (ie,
	% k == 2). In that case, it goes from some lower bound (jIII) up to nK.
	%
	jmesh = jIII:nK;
	
	%
	% Loop over the two Phases. Phase 1 = go up in K2. Phase 2: start at
	% some K2 and go down to the lowest value of K2.
	%
	flipped = false;
	for Phase = 1:2
		if Phase == 1
			imesh = 1:nK;
		else
			imesh = iII:-1:1;
		end
		
		%
		% Loop over K2.
		%
		for ii = 1:length(imesh)
			i = imesh(ii);
			K2 = K22(i);
			
			%
			% At some point in incrementing up K2, we decide to switch the
			% direction that we loop over K1 (the inner loop), so jmesh
			% changes.
			%
			if ~flipped && K2 > 1e-2%K2_I_FBL(1)
				jmesh = nK:-1:jIII;
				xp0 = log(1.0620e-05);
% 				xp0 = log(0.0615);
				flipped = true; % "flipped" variable only here so that this
				% code to flip jmesh only gets executed once. This is so
				% that the initial guess on x only gets applied once.
			end
			
			%
			% Initial guess is either something totally new, or the
			% prediction, xp0, from the last value of i (or, one time only,
			% it's the flip).
			%
			if i == imesh(1)
				if k == 1
					x = log(K12_PA(imesh(1),jmesh(1)));
				else
					x = 0;
				end
			else
				x = xp0;
			end
			
			%
			% Define the stopping criterion: if K2 < max(K2_II_FBL), then
			% we will run into Region II as we move in K1, so that will be
			% our break-out criterion.
			%
			if i <= iII && Phase == 1
				K1max = spline(K2_II_FBL(iII_max:end),K1_II_FBL(iII_max:end),K2);
			elseif i <= iII && Phase == 2
				K1min = spline(K2_II_FBL(1:iII_max),K1_II_FBL(1:iII_max),K2);
			elseif i >= iI && Phase == 1
				K1min = spline(K2_I_FBL,K1_I_FBL,K2);
			end
			
			%
			% Continuing along the K1 direction
			%
			for jj = 1:length(jmesh)
				j = jmesh(jj);
				% Stopping criterion: if we hit Region II.
				K1 = K11(j); K12 = K12_PA(i,j);
				if (i <= iII && ((Phase == 1 && K1 > K1max) || (Phase == 2 && K1 < K1min))) ...
						|| (i > iI && Phase == 1 && K1 < K1min)
					break
				end
				p = [1 K1 K2 K12 K3 K4 1 1];				
				
				%
				% While loop for Newton's method to find the value of K12
				% that give NPA.
				%
				f = 1; dx = 0; nSteps = 0;
				tolerance = 1e-6; nStepsMax = 25;
				while norm([f;dx]) > tolerance && nSteps < nStepsMax
					
					%
					% calc f
					%
					p(4) = exp(x);
					[~,~,w0] = ftn_goentoro_ss(p,x0,n);
					[~,~,w1] = ftn_goentoro_ss(p,x1,n);
					A = (1 + w0/K3)./(1 + w1/K3);
					f = fhandle(K1,K2,exp(x),A);
					
					%
					% Calc fprime, the derivative of f wrt K12
					%
					x_ = x + delt;
					p1 = p;
					p1(4) = exp(x_);
					[~,~,w0] = ftn_goentoro_ss(p1,x0,n);
					[~,~,w1] = ftn_goentoro_ss(p1,x1,n);
					A = (1 + w0/K3)./(1 + w1/K3);
					f1 = fhandle(K1,K2,exp(x_),A);
					dfdlogK12 = (f1 - f)/delt;
					
					%
					% Take the step
					%
					dx = -f/dfdlogK12;
					x = x + dx;
					nSteps = nSteps + 1;
				end
				K12_NPA_FBL(i,j,k) = exp(x); % output
				
				%
				% Saving the predicted x for when we increment K2.
				%
				if j == jmesh(1) && ii < length(imesh)
					
					%
					% Now calc. the derivative of f wrt K2
					%
					K2_1 = K2*(1 + delt);
					p1 = p;
					p1(3) = K2_1;
					[~,~,w0] = ftn_goentoro_ss(p1,x0,n);
					[~,~,w1] = ftn_goentoro_ss(p1,x1,n);
					A = (1 + w0/K3)./(1 + w1/K3);
					f1 = fhandle(K1,K2_1,exp(x),A);
					dfdlogK2 = (f1 - f)/delt;
					dlogK12dlogK2 = - dfdlogK2/dfdlogK12;
					dlogK2 = log(K22(imesh(ii+1))) - log(K2);
					xp0 = x + dlogK12dlogK2*dlogK2; % next guess
				end
				
				%
				% Next predicted x for when we increment K1
				%
				if jj < length(jmesh)
					
					%
					% Now calc. the derivative of f wrt K1
					%
					K1_1 = K1*(1 + delt);
					p1 = p;
					p1(2) = K1_1;
					[~,~,w0] = ftn_goentoro_ss(p1,x0,n);
					[~,~,w1] = ftn_goentoro_ss(p1,x1,n);
					A = (1 + w0/K3)./(1 + w1/K3);
					f1 = fhandle(K1_1,K2,exp(x),A);
					dfdlogK1 = (f1 - f)/delt;
					dlogK12dlogK1 = - dfdlogK1/dfdlogK12;
					dlogK1 = log(K11(jmesh(jj+1))) - log(K1);
					x = x + dlogK12dlogK1*dlogK1; % next guess
				end
				
			end
% 			disp(['i = ',num2str(i),' out of ',num2str(nK),'; k = ',num2str(k)])
		end
		disp(['Phase = ',num2str(Phase),' of 2; k = ',num2str(k)])
	end
end

save Mat/Fig3_FBL_K12NPA_K1K2sweep




