% script_FigSXX_w_coop
%
% This script is designed to run all simulations, etc, to get Fig 3 of
% Reeves 2019 (FF/FB).

clear
close all
options = ddeset('RelTol',1e-6);

%
% Param values
%
x0 = 1;
x1 = 10;
n = 1;
vareps = 0.05;
delt = 1e-4;

thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
tauz = 1; tauw = 1; tauy = 1;
K3 = 0.1; K4 = 0.1;

%
% Disturbance function
%
t0 = 0;
tspan = [t0 100];
disthand = @(t)ftn_unitstep(t,t0,x0,x1);

CC13 = [1 0.5 0.3 0.1];
CC23 = [1 0.5 0.3 0.1];




%% ========================================================================
% First, we will calculate the Delta(K12)/K12 for the situation with
% negative feedback. We are also adding the dde simulations. We have to combine
% them because we are also varying C13,C23, so all of our plots, both the
% K12_NPA AND the dde plots, will have to be with these varied C13,C23.
% =========================================================================
%{

% load Mat/DeltaK_FFLonly_dense K11_dense K22_dense DeltaK
% load Mat/C_pointohone_contour_coop XC YC
% load Mat/Peaks_pointone_contours_FFFB XP YP Xzmax Yzmax VP Vzmax
% load Mat/Peaks_pointone_contours_FFFB_w_coop XP YP Xzmax Yzmax VP Vzmax


for i13 = 1:length(CC13)
	C13 = CC13(i13);
	for i23 = 1:length(CC23)
		C23 = CC23(i23);
		
		
		
		% -----------------------------------------------------------------
		% First, the NPA calculations
		% -----------------------------------------------------------------
		
		%
		% Param values
		%		
		nK = 500;
		K11 = logspace(-2,6,nK)';
		K22 = logspace(-6,2,nK)';
		
	
		%
		% creating K12_PA and the K12_NPAplus and minus arrays
		%
		script_regionsI_II_coop
		script_FigS4_K12NPA_K1K2sweep
		
		K11_NPA = K11; K22_NPA = K22;


		
		% -----------------------------------------------------------------
		% Now, the dde calculations
		% -----------------------------------------------------------------
		
		%
		% Perfect adaptation, on the less-dense mesh
		%
		nK = 50;
		K11 = logspace(-2,6,nK)';
		K22 = logspace(-6,2,nK)';
		y0 = x0./(K11 + x0);
		y1 = x1./(K11 + x1);
		RHS = y1' - y0';
		LHS = (1 + (1./K22)*y0')/x0 - (1 + (1./K22)*y1')/x1;
		K12_PA = repmat(RHS,nK,1)./LHS;		
		
		%
		% All (slow) simulation commands in the following script, which is
		% commented out now, since the sims were already run and stored in the Mat
		% file below.
		%
% 		script_FigS4_dde_K1K2sweep
		
		data.x0 = x0;
		data.x1 = x1;
		data.n = n;
		data.vareps = vareps;
		data.K3 = K3;
		data.K4 = K4;
		data.C13 = C13;
		data.C23 = C23;

		data.theta = theta;
		data.tau = [tauz tauw tauy];
		
% 		data.K1_span = {log10(K11(1)) log10(K11(end)) length(K11)};
% 		data.K2_span = {log10(K22(1)) log10(K22(end)) length(K22)};
% 		data.Zpeak = Zpeak;
% 		data.Zinitial = Zinitial;
		
		
		data.K1_NPA_span = {log10(K11_NPA(1)) log10(K11_NPA(end)) length(K11_NPA)};
		data.K2_NPA_span = {log10(K22_NPA(1)) log10(K22_NPA(end)) length(K22_NPA)};
		data.K12_NPA_FBL = K12_NPA_FBL;
		data.K12_PA_FBL = K12_PA_FBL;
		
		if i13 == 1 && i23 == 1
			fnames = fieldnames(data);
			CELL_pre = cell(length(CC13),length(CC23));
			CELL_pre2 = [fnames repmat({CELL_pre},length(fnames),1)]';
			Soln2 = struct(CELL_pre2{:});
		end
		Soln2(i13,i23) = data;

		disp(['i13 = ',num2str(i13),' out of ',num2str(length(CC13)),...
			', i23 = ',num2str(i23),' out of ',num2str(length(CC23))])
	end
end
save Mat/FigS4_Soln_all2 Soln2

%}


%% ========================================================================
% Now that we have run simulations, and also found K12_PA and NPA, we can
% plot heatmaps.
% =========================================================================
% {


load Mat/FigS4_Soln_all Soln
load Mat/FigS4_Soln_all2 Soln2

%
% Load FF only contours (C_PA = 0.01, P = 0.1, Zmax = 0.01)
%
load Mat/C_pointohone_contour_coop XC YC
% load Mat/Peaks_pointone_contours_FFFB XP YP Xzmax Yzmax
load Mat/Peaks_pointone_contours_FFFB X_realistic Y_realistic

% -------------------------------------------------------------------------
% Calculate the NPA for FF only
% -------------------------------------------------------------------------

nK = Soln2(1,1).K1_NPA_span{3};
K11 = logspace(Soln2(1,1).K1_NPA_span{:})';
K22 = logspace(Soln2(1,1).K2_NPA_span{:})';
K111 = repmat(K11',nK,1);

%
% K12_PA
%
x0 = 1;
x1 = 10;
y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
RHS = y1' - y0';
LHS = (1/x0-1/x1) - (1./K22)*(1./(K11' + x1) - 1./(K11' + x0));
K12_PA = repmat(RHS,nK,1)./LHS;

%
% creating K12_NPAplus and minus arrays
%
RHS = K11'.*(y1' - y0'/(1 + vareps));
LHS = 1/(1 + vareps)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAplus = repmat(RHS,nK,1)./LHS;

% K12_NPAminus
vareps1 = -vareps;
RHS = K11'.*(y1' - y0'/(1 + vareps1));
LHS = 1/(1 + vareps1)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAminus = repmat(RHS,nK,1)./LHS;

K1lower_minus = vareps1/(1/x1 - 1/x0*(1 + vareps1));
K12_NPAminus(K11 <= K1lower_minus) = 0;
K12_NPAminus(K12_NPAminus < 0) = 0;

Y = (K12_NPAplus-K12_PA)./K12_PA;
Yminus = (K12_PA-K12_NPAminus)./K12_PA;
DeltaK = Y + Yminus;



% -----------------------------------------------------------------
% Plotting of the dde stuff, get the contours
% -----------------------------------------------------------------

figure('position',[337.0000   61.8000  946.4000  707.2000])
set(gcf,'paperpositionmode','auto')
for i13 = 1:length(CC13)
	C13 = CC13(i13);
	for i23 = 1:length(CC23)
		C23 = CC23(i23);

		data = Soln(i13,i23);
		data2 = Soln2(i13,i23);
		
		% -----------------------------------------------------------------
		% Get the contours that come from the dde stuff
		% -----------------------------------------------------------------
		
		%
		% Extract the dde data from the structure.
		%
		nK = data.K1_span{3};
		K11 = logspace(data.K1_span{:})';
		K22 = logspace(data.K2_span{:})';
		Zpeak = data.Zpeak;
		Zinitial = data.Zinitial;
		
		%
		% Make super-dense mesh
		%
		nK2 = 1000;
		K11_dense2 = logspace(data.K1_span{1:2},nK2)';
		K22_dense2 = logspace(data.K2_span{1:2},nK2)';
		
		%
		% Getting the 0.1 contour for P
		%
		P = (Zpeak - Zinitial)./Zinitial;
% 		c = contourc(K11,K22,P,[0.1 0.1]);
% 		[X1,Y1] = extractcontour(c);
% 		N = cellfun(@length,X1);
% 		[~,imax] = max(N);
% 		XP_coop = X1{imax};
% 		YP_coop = Y1{imax};
		P2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),P,...
			repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
		VP = P2 >= 0.1;
				
		%
		% Getting the 0.01 contour for Zpeak
		%		
% 		c = contourc(K11,K22,Zpeak,[0.01 0.01]);
% 		[X1,Y1] = extractcontour(c);
% 		N = cellfun(@length,X1);
% 		[~,imax] = max(N);
% 		Xzmax_coop = X1{imax};
% 		Yzmax_coop = Y1{imax};
		Zpeak2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),Zpeak,...
			repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
		Vzmax = Zpeak2 >= 0.01;
		
		% -----------------------------------------------------------------
		% Finally, the plotting for the NPA stuff
		% -----------------------------------------------------------------
		
		%
		% Extract the NPA data from the structure.
		%
		nK = data2.K1_NPA_span{3};
		K11 = logspace(data2.K1_NPA_span{:})';
		K22 = logspace(data2.K2_NPA_span{:})';
		K12_PA_FBL = data2.K12_PA_FBL;
		K12_NPA_FBL = data2.K12_NPA_FBL;
		script_regionsI_II_coop
		
		
		%
		% Getting the 0.01 contour for C_PA
		%
		C_coop = K12_PA_FBL./repmat(K11',nK,1)./repmat(K22,1,nK);
% 		c = contourc(K11,K22,C_coop,[0.01 0.01]);
% 		[X1,Y1] = extractcontour(c);
% 		N = cellfun(@length,X1);
% 		[~,imax] = max(N);
% 		XC_coop = X1{imax};
% 		YC_coop = Y1{imax};
		C_coop2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),C_coop,...
			repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
		VC = C_coop2 >= 0.01;
		
		%
		% Realistic region for coop
		%
		V_realistic = VP & VC & Vzmax;
		if any(V_realistic(:))
			[~,J_realistic,I_realistic] = traceobject(V_realistic);
			X_realistic_coop = K11_dense2(J_realistic);
			Y_realistic_coop = K22_dense2(I_realistic);
		else
			X_realistic_coop = [];
			Y_realistic_coop = [];
			
		end
		
		
		
		%
		% Plotting DeltaK12^{tot}/K12
		%
		Y = (real(K12_NPA_FBL(:,:,1)) - K12_PA_FBL)./K12_PA_FBL;
		Yminus = (K12_PA_FBL - real(K12_NPA_FBL(:,:,2)))./K12_PA_FBL;
		DeltaK_FBL = Y + Yminus;
% 		figure('paperpositionmode','auto')
% 		pcolor_contour(K11,K22,DeltaK_FBL,[0.15 0.3 1 5],'r',true,0,5)
% 		set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
% 		
% 		plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
% 		plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
% 		plot(XC_coop,YC_coop,'w',XP_coop,YP_coop,'w',Xzmax_coop,Yzmax_coop,'w')
% 		plot(XC,YC,'m',XP,YP,'m',Xzmax,Yzmax,'m')
% 		savepcolor(gcf,'Figs/DeltaK12_NPA_FFFB_w_coop')
% 		xlabel('K1')
% 		ylabel('K2')
% 		title('\Delta K_3^{NPA,FBL}/K_3^{PA}')
		
		
		

		
		%
		% pcolor and contour of the ratio of the two DeltaK's. Note that in this
		% contour, there is a weird curve for the ratio equal to 1.5.
		%
		RDK = DeltaK_FBL./DeltaK; % ratio of DeltaK's
% 		figure('paperpositionmode','auto')
		subplot(length(CC13),length(CC23),(i13-1)*length(CC23)+i23)
		pcolor_contour(K11,K22,RDK,[1 1.1 1.25 1.5 1.75 2],'r',true,0.99,2,[])
		set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
		
		plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
		plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
% 		plot(XC_coop,YC_coop,'w',XP_coop,YP_coop,'w',Xzmax_coop,Yzmax_coop,'w')
% 		plot(XC,YC,'m',XP,YP,'m',Xzmax,Yzmax,'m')
		plot(X_realistic_coop,Y_realistic_coop,'w')	
		plot(X_realistic,Y_realistic,'k')	
% 		
% 		h = get(gcf,'children');
% 		for i = 1:length(h)
% 			if strcmp(h(i).Tag,'Colorbar')
% 				cbh = h(i);
% 				set(h(i),'YTick',1:0.25:2)
% 			end
% 		end
% 		savepcolor(gcf,'Figs/Ratio_of_DeltaK12_w_coop')
		if i13 == length(CC13)
			xlabel('K1')
		end
		if i23 == 1
			ylabel('K2')
		end
		title(['C_{13} = ',num2str(C13),', C_{23} = ',num2str(C23)])

	end
end


h = get(gcf,'children');
% tobedeleted = false(length(h),1);
for i = 1:length(h)
	if strcmp(h(i).Tag,'Colorbar')
		cbh = h(i);
		set(h(i),'YTick',1:0.25:2)
	end
		
	if strcmp(h(i).Tag,'legend')
% 		tobedeleted(i) = true;
		delete(h(i))
	end
end
% for i = 1:length(h)
% 	if tobedeleted(i)
% 		delete(h(i))
% 	end
% end

print(gcf,'Figs/FigS4_coop.jpg','-djpeg','-r300')


%}
