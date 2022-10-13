% script_FigS4_K12NPA_K1K2sweep
%
% This script goes through each of the values of K1,K2 and finds the value
% of K3 that gives NPA. We use a double-for loop and a Newton in the
% middle. We will choose moderate levels of neg fbk strength.



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
K12_NPA_FBL = NaN(nK,nK,3);

%
% Need max of and index of max of K2 along Region II.
%
[K2_II_FBL_max,iII_max] = max(K2_II_FBL);
iII = find(K22 < K2_II_FBL_max);
iII = iII(end);

%
% Need index of min of K2 along Region I.
%
iI = find(K22 <= K2_I_FBL(1));
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
for k = 1:3
	if k == 2
		vareps1 = -vareps1;
		
		%
		% If vareps < 0, then there is a lower limit to K1, below which K12
		% can go all the way down to zero and you still get NPA. IOW, it
		% doesn't make sense to find the lower bound on K12 for NPA if K1
		% is too low. So we'll skip that part.
		%
		K1_lower = vareps1/(1/x1 - 1/x0*(1 + vareps1));
		
		%
		% Establish this K1_lower as a lower bound.
		%
		jIII = find(K11 <= K1_lower);
		jIII = jIII(end) + 1;
		
		K12_NPA_FBL(:,1:jIII-1,k) = 0;
		
	elseif k == 3
		vareps1 = 0;
		jIII = 1;	
	else
		jIII = 1;	
	end
	
	%
	% Define function handles.
	%
	fhandle = @ftn_FB_K1K2_sweep_w_coop;
	
	%
	% Loop over the two Phases. Phase 1 = go up in K2. Phase 2: start at
	% some K2 and go down to the lowest value of K2.
	%
	% Here we also define jmesh, for varying K1. This goes from nK down to
	% 1, unless you have vareps < 0 (ie, k == 2). In that case, it goes
	% from nK down to some lower bound (jIII).
	%
	flipped = false;
	for Phase = 1:2
		if Phase == 1
			imesh = 1:nK;
			jmesh = nK:-1:jIII;
		else
			imesh = iII:-1:1;
			jmesh = jIII:nK;
		end
		
		%
		% Loop over K2.
		%
		for ii = 1:length(imesh)
			i = imesh(ii);
			K2 = K22(i);
			
			%
			% Initial guess is either something totally new, or the
			% prediction, xp0, from the last value of i (or, one time only,
			% it's the flip).
			%
			if i == imesh(1)
				if Phase == 2
					x = xp01;
				else
					x = log(1e-5);
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
				K1min = pchip(K2_II_FBL(1:iII_max),K1_II_FBL(1:iII_max),K2);
% 				K1min = spline(K2_II_FBL(1:iII_max),K1_II_FBL(1:iII_max),K2);
			elseif i <= iII && Phase == 2
				K1max = pchip(K2_II_FBL(iII_max:end),K1_II_FBL(iII_max:end),K2);
% 				K1max = spline(K2_II_FBL(iII_max:end),K1_II_FBL(iII_max:end),K2);
			elseif i >= iI && Phase == 1
				K1min = pchip(K2_I_FBL,K1_I_FBL,K2);
% 				K1min = spline(K2_I_FBL,K1_I_FBL,K2);
			end
			
			%
			% Continuing along the K1 direction
			%
			for jj = 1:length(jmesh)
				j = jmesh(jj);
				% Stopping criterion: if we hit Region II.
				K1 = K11(j); %K12 = K12_PA(i,j);
				if (i <= iII && ((Phase == 1 && K1 < K1min) || ...
						(Phase == 2 && K1 > K1max))) ...
						|| (i > iI && Phase == 1 && K1 < K1min)
					break
				end
% 				p = [1 K1 K2 K12 K3 K4 1 1 C13 C23];	
				p = [1 K1 K2 1 K3 K4 1 1 C13 C23];				
				
				%
				% While loop for Newton's method to find the value of K12
				% that give NPA.
				%
				x_0 = x;
				f = 1; dx = 0; nSteps = 0;
				tolerance = 1e-6; nStepsMax = 25;
				while norm([f;dx]) > tolerance && nSteps < nStepsMax
					
					%
					% calc f
					%
					p(4) = exp(x);
					f = fhandle(p,x0,x1,n,vareps1);
					
					%
					% Calc fprime, the derivative of f wrt K12
					%
					x_1 = x + delt;
					p1 = p;
					p1(4) = exp(x_1);
					f1 = fhandle(p1,x0,x1,n,vareps1);
					dfdlogK12 = (f1 - f)/delt;
					
					%
					% Take the step
					%
					dx = -f/dfdlogK12;
					x = x + dx;
					nSteps = nSteps + 1;
				end
				
				%
				% Checking the output of Newton
				%
				if exp(x) < 1e-16
					K12_NPA_FBL(i,jmesh(jj:end),k) = 0;
					break
					
					
				elseif isnan(x) || nSteps == nStepsMax % added the OR
					%
					% If our Newton did not converge, we will try to
					% converge by brute-force + advNtn.
					%
					fhandle1 = @(x)fhandle([p(1:3) exp(x) p(5:end)],x0,x1,n,vareps1);
					KD = logspace(x_0/log(10)-10,x_0/log(10)+5,200)';
					FKD = zeros(length(KD),1);
					for ikd = 1:length(KD)
						FKD(ikd) = fhandle1(log(KD(ikd)));
					end
					
					k0 = find((FKD(1:end-1) < 0 & FKD(2:end) > 0) |...
						(FKD(1:end-1) > 0 & FKD(2:end) < 0));
					[~,imin] = min(abs(KD(k0)-(x_0/log(10))));
					k0 = k0(imin);
					
					if ~isempty(k0)
						[x,dfdlogK12,f] = advNtnDU(fhandle1,x_0,...
							[log(KD(k0:k0+1))';FKD(k0:k0+1)']); % added "f" output
						if ischar(x)
							K12_NPA_FBL(i,jmesh(jj:end),k) = 0;
							break
						else
							
							K12_NPA_FBL(i,j,k) = exp(x); % output
							p(4) = exp(x); % added update to "p"
						end
					elseif Phase == 1
						K12_NPA_FBL(i,jmesh(jj:end),k) = 0;
						break
					elseif Phase == 2
						break
					end					
					
					
					
				else
					K12_NPA_FBL(i,j,k) = exp(x); % output
					p(4) = exp(x); % added update to "p"
				end
				
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
					f1 = fhandle(p1,x0,x1,n,vareps1);
					
					dfdlogK2 = (f1 - f)/delt;
					dlogK12dlogK2 = - dfdlogK2/dfdlogK12;
					dlogK2 = log(K22(imesh(ii+1))) - log(K2);
					xp0 = x + dlogK12dlogK2*dlogK2; % next guess
				
				%
				% Saving a predicted x for when we flip and do Phase 2
				%
				elseif Phase == 1 && i == iII+1 && j == jmesh(end)
					%
					% Now calc. the derivative of f wrt K2
					%
					K2_1 = K2*(1 + delt);
					p1 = p;
					p1(3) = K2_1;
					f1 = fhandle(p1,x0,x1,n,vareps1);
					
					dfdlogK2 = (f1 - f)/delt;
					dlogK12dlogK2 = - dfdlogK2/dfdlogK12;
					dlogK2 = log(K22(imesh(ii+1))) - log(K2);
					xp01 = x - dlogK12dlogK2*dlogK2; % next guess
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
					f1 = fhandle(p1,x0,x1,n,vareps1);
					dfdlogK1 = (f1 - f)/delt;
					dlogK12dlogK1 = - dfdlogK1/dfdlogK12;
					dlogK1 = log(K11(jmesh(jj+1))) - log(K1);
					x = x + dlogK12dlogK1*dlogK1; % next guess
				end
				
			end
		end
		disp(['Completed Phase ',num2str(Phase),' of 2; k = ',num2str(k)])
	end
end

K12_PA_FBL = K12_NPA_FBL(:,:,3);
K12_NPA_FBL(:,:,3) = [];
C = K12_PA_FBL./repmat(K11',nK,1)./repmat(K22,1,nK);


% save Mat/FigS4_FBL_K12NPA_K1K2sweep_w_coop

%%

%{

%
% Plotting, for the heck of it
%
veps = vareps*[1 -1];
for k = 1:2
	figure
	pcolor_contour(K11,K22,log10(K12_NPA_FBL(:,:,k)),[],'',true,-3,1)
	loglog(K1_I_FBL,K2_I_FBL)
	loglog(K1_II_FBL,K2_II_FBL)
	title(['\epsilon = ',num2str(veps(k)),', C_{13} = ',...
		num2str(C13),', C_{23} = ',num2str(C23)])
end
% 
% figure
% pcolor_contour(K11,K22,log10(K12_PA),[],'',true,-3,1)
% loglog(K1_I_FBL,K2_I_FBL)
% loglog(K1_II_FBL,K2_II_FBL)
% title(['\epsilon = ',num2str(0),', C_{13} = ',num2str(C13),', C_{23} = ',num2str(C23)])



% %
% % Getting the 0.01 contour for C_PA
% %
% c = contourc(K11,K22,C,[0.01 0.01]);
% [X1,Y1] = extractcontour(c);
% for i = 1:length(X1)
% 	if length(X1{i}) > 100
% 		XC = X1{i};
% 		YC = Y1{i};
% 		break
% 	end
% end
% 
% % figure
% % pcolor_contour(K11,K22,C,[0.01 0.1 0.5 1 1.5 2 2.5],'',true)
% % loglog(K1_I_FBL,K2_I_FBL)
% % loglog(K1_II_FBL,K2_II_FBL)
% % title(['\epsilon = 0, C_{13} = ',num2str(C13),', C_{23} = ',num2str(C23)])
% % 
% % %
% % % Get the C = 0.01 contour
% % %
% % h = get(gca,'children');
% % for i = 1:length(h)
% % 	if strcmp(h(i).DisplayName,'0.01')
% % 		XC = get(h(i),'Xdata')';
% % 		YC = get(h(i),'Ydata')';
% % 		break
% % 	end
% % end
% 
% save Mat/C_pointohone_contour_coop XC YC C13 C23 K3 K4


%}



