% script_FigS6_FFL_dde_K1K2sweep



%
% Make a for loop to run simulations for many different K1,K2's
%
Zinitial_2 = zeros(nK,nK);
Zpeak_2 = zeros(nK,nK);
Tpeak_2 = zeros(nK,nK);
Zfinal_2 = zeros(nK,nK);
Zfinal_21 = zeros(nK,nK);
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
		y0 = x0_2./(K1 + x0_2);
		z0 = x0_2/K1./(1 + x0_2/K1 + y0/K2 + x0_2*y0/K12);
		w0 = 0;
		
		%
		% Final conditions
		%
		y1 = x1_2./(K1 + x1_2);
		z1 = x1_2/K1./(1 + x1_2/K1 + y1/K2 + x1_2*y1/K12);
		
		
		%
		% FF only sim
		%
		ftnhand = @ftn_goentoro_dde;
		Y0 = [y0; z0; w0];
		soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand2);
% 		soln = ftn_rundde(p,theta,x0_2,x1_2,n);
		t = soln.x';
		Y = soln.y';
		z = Y(:,2);
		Zinitial_2(i,j) = z0;
		[Zpeak_2(i,j),ipeak] = max(z);
		Tpeak_2(i,j) = t(ipeak) - t0;
		Zfinal_2(i,j) = z1;
		
		
	end
	disp(['i = ',num2str(i),' out of ',num2str(nK)])
end
save Mat/FigS6_peak_K1K2sweep