% script_FigS4_dde_K1K2sweep


% load Mat/FigSXX_FBL_K12NPA_K1K2sweep_w_coop K12_PA C13 C23 K3 K4


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
		p = [tauz K1 K2 K12 K3 K4 tauw tauy C13 C23];
		
		if K1 > 1.7e4
			
			options = ddeset('RelTol',1e-8,'AbsTol',1e-8);
		else
			
			options = ddeset('RelTol',1e-6);
		end
		
		% =============================================================
		% FFFB
		% =============================================================
		
		%
		% Initial conditions
		%
		[y0,z0,w0] = ftn_goentoro_ss(p,x0,n);
		
		%
		% FF/FB sim
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
		Zfinal(i,j) = z(end);
		
		
	end
	disp(['i = ',num2str(i),' out of ',num2str(nK)])
end
save Mat/FigSXX_peak_K1K2sweep_FFFB_coop