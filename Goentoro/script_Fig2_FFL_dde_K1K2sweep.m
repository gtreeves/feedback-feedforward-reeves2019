% script_Fig2_FFL_dde_K12sweep

y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
RHS = y1' - y0';
LHS = (1 + (1./K22)*y0')/x0 - (1 + (1./K22)*y1')/x1;
K12_PA = repmat(RHS,nK,1)./LHS;

%
% Make a for loop to run simulations for many different K1,K2's
%
Zinitial = zeros(nK,nK);
Zpeak = zeros(nK,nK);
Tpeak = zeros(nK,nK);
Zfinal = zeros(nK,nK);
for i = 1:nK
	K2 = K22(i); 
	
	for j = 1:nK
		K1 = K11(j); K12 = K12_PA(i,j);
		K3 = 1; K4 = Inf;
		p = [tauz K1 K2 K12 K3 K4 tauw tauy];
		
		if K1 > 1.7e4
			
			options = ddeset('RelTol',1e-8,'AbsTol',1e-8);
		else
			
			options = ddeset('RelTol',1e-6);
		end
		
		% =============================================================
		% FF only
		% =============================================================
		
		%
		% Initial conditions
		%
		y0 = x0./(K1 + x0);
		z0 = x0/K1./(1 + x0/K1 + y0/K2 + x0*y0/K12);
		w0 = 0;
		
		%
		% Final conditions
		%
		y1 = x1./(K1 + x1);
		z1 = x1/K1./(1 + x1/K1 + y1/K2 + x1*y1/K12);
		
		%
		% FF only sim
		%
		ftnhand = @ftn_goentoro_dde;
		Y0 = [y0; z0; w0];
		soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
		t = soln.x';
		Y = soln.y';
		z = Y(:,2);
		Zinitial(i,j) = z0;
		[Zpeak(i,j),ipeak] = max(z);
		Tpeak(i,j) = t(ipeak) - t0;
		Zfinal(i,j) = z1;
		
		
	end
	disp(['i = ',num2str(i),' out of ',num2str(nK)])
end
save Mat/Fig2_peak_K1K2sweep