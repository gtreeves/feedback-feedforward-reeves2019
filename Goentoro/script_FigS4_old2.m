% script_FigSXX_w_coop
%
% This script is "old2" because we used to plot all of the dynamic stuff,
% but now we're just doing contourc to get the contours.
%
% This script is designed to run all simulations, etc, to get Fig 3 of
% Reeves 2019 (FF/FB).

% clear
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


for i13 = 2:length(CC13)
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
		script_FigS4_dde_K1K2sweep
		
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
		
		data.K1_span = {log10(K11(1)) log10(K11(end)) length(K11)};
		data.K2_span = {log10(K22(1)) log10(K22(end)) length(K22)};
		data.Zpeak = Zpeak;
		data.Zinitial = Zinitial;
		
		
		data.K1_NPA_span = {log10(K11_NPA(1)) log10(K11_NPA(end)) length(K11_NPA)};
		data.K2_NPA_span = {log10(K22_NPA(1)) log10(K22_NPA(end)) length(K22_NPA)};
		data.K12_NPA_FBL = K12_NPA_FBL;
		data.K12_PA_FBL = K12_PA_FBL;
		
		if i13 == 1 && i23 == 1
			fnames = fieldnames(data);
			CELL_pre = cell(length(CC13),length(CC23));
			CELL_pre2 = [fnames repmat({CELL_pre},length(fnames),1)]';
			Soln = struct(CELL_pre2{:});
		end
		Soln(i13,i23) = data;

		disp(['i13 = ',num2str(i13),' out of ',num2str(length(CC13)),...
			', i23 = ',num2str(i23),' out of ',num2str(length(CC23))])
	end
end
save Mat/FigS4_Soln_all Soln

%}


%% ========================================================================
% Now that we have run simulations, and also found K12_PA and NPA, we can
% plot heatmaps.
% =========================================================================
% {


load Mat/FigS4_Soln_all Soln

%
% Contour of C_PA = 0.01
%
load Mat/C_pointohone_contour_coop XC YC



% -----------------------------------------------------------------
% Plotting of the dde stuff, get the contours
% -----------------------------------------------------------------

for i13 = 1:length(CC13)
	C13 = CC13(i13);
	for i23 = 1:length(CC23)
		C23 = CC23(i23);

		data = Soln(i13,i23);
		Zpeak = data.Zpeak;
		Zinitial = data.Zinitial;

		%
		% Plotting pcolor and contour of P.
		%
		P = (Zpeak - Zinitial)./Zinitial;
		figure('paperpositionmode','auto')
		pcolor_contour(K11,K22,P,[0.01 0.1 1 1.95 2.05 3],'r',true)
		set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
		plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
		plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
		% savepcolor(gcf,'Figs/Peak_surface_FFFB')
		xlabel('K1')
		ylabel('K2')
		title('P = (Z_{max}-Z_0)/Z_0')
		
		%
		% Get the 0.1 contour
		%
		h = get(gca,'children');
		for i = 1:length(h)
			if strcmp(h(i).DisplayName,'0.1')
				XP = get(h(i),'Xdata')';
				YP = get(h(i),'Ydata')';
				break
			end
		end

		%
		% Plotting pcolor and contour of Zpeak. This is for Fig. S2.
		%
		figure('paperpositionmode','auto')
		pcolor_contour(K11,K22,Zpeak,[0.01 0.1 1],'r',true)
		set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
		plot(XP,YP,'w')
		plot(XC,YC,'w')
		plot(K1_I_FBL,K2_I_FBL,'k','linewidth',1)
		plot(K1_II_FBL,K2_II_FBL,'k','linewidth',1)
		% plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
		% savepcolor(gcf,'Figs/Zmax_surface_FFFB')
		xlabel('K1')
		ylabel('K2')
		title('Z_{max}')
		
		%
		% Get the zmax = 0.01 contour
		%
		h = get(gca,'children');
		for i = 1:length(h)
			if strcmp(h(i).DisplayName,'0.01')
				Xzmax = get(h(i),'Xdata')';
				Yzmax = get(h(i),'Ydata')';
				break
			else
				Xzmax = [];
				Yzmax = [];
			end
		end
		
		
		
		% -----------------------------------------------------------------
		% Finally, the plotting for the NPA stuff
		% -----------------------------------------------------------------
		
		
		
		%
		% Plotting DeltaK12^{tot}/K12
		%
		Y = (real(K12_NPA_FBL(:,:,1)) - K12_PA)./K12_PA;
		Yminus = (K12_PA - real(K12_NPA_FBL(:,:,2)))./K12_PA;
		DeltaK_FBL = Y + Yminus;
		figure('paperpositionmode','auto')
		pcolor_contour(K11,K22,DeltaK_FBL,[0.15 0.3 1 5],'r',true,0,5)
		set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
		
		plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
		plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
		plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
		% plot(XC,YC,'w',XP,YP,'w',Xzmax,Yzmax,'w')
		savepcolor(gcf,'Figs/DeltaK12_NPA_FFFB_w_coop')
		xlabel('K1')
		ylabel('K2')
		title('\Delta K_3^{NPA,FBL}/K_3^{PA}')
		
		
		DeltaK1 = interp2(repmat(K11_dense',length(K22_dense),1),...
			repmat(K22_dense,1,length(K11_dense)),DeltaK,...
			repmat(K11',nK,1),repmat(K22,1,nK));
		
		%
		% pcolor and contour of the ratio of the two DeltaK's. Note that in this
		% contour, there is a weird curve for the ratio equal to 1.5.
		%
		RDK = DeltaK_FBL./DeltaK1; % ratio of DeltaK's
		figure('paperpositionmode','auto')
		pcolor_contour(K11,K22,RDK,[1 1.1 1.25 1.5 1.75 2],'r',true,0.99,2)
		set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
		plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
		plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
		
		plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
		% plot(XC,YC,'w',XP,YP,'w',Xzmax,Yzmax,'w')
		h = get(gcf,'children');
		for i = 1:length(h)
			if strcmp(h(i).Tag,'Colorbar')
				cbh = h(i);
				set(h(i),'YTick',1:0.25:2)
			end
		end
		savepcolor(gcf,'Figs/Ratio_of_DeltaK12_w_coop')
		xlabel('K1')
		ylabel('K2')
		title('\Delta K_3^{NPA,FBL}/\Delta K_3^{NPA}')

	end
end

%}



%% ========================================================================
% Now that we have found the values of K12 for which f  = eps = 0.05 (and
% for f = -0.05, we will run a set of simulations for a particular choice
% of K1 and K2 to show what it looks like.
% =========================================================================
%{

clear
options = ddeset('RelTol',1e-6);
load Mat/Fig3_FBL_K12NPA_K12sweep
K1 = 1; K2 = 0.1; 

K12_NPAplus = interp2(repmat(K11',nK,1),repmat(K22,1,nK),K12_NPA_FBL(:,:,1),K1,K2);
K12_NPAminus = interp2(repmat(K11',nK,1),repmat(K22,1,nK),K12_NPA_FBL(:,:,2),K1,K2);

%
% K12_PA
%
F0 = 1;
F1 = 10;
y0 = F0./(K1 + F0);
y1 = F1./(K1 + F1);
RHS = y1 - y0;
LHS = (1/F0-1/F1) - (1./K2)*(1./(K1 + F1) - 1./(K1 + F0));
K12 = RHS/LHS;



%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Param values
%
n = 1;
thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 20];
K3 = 0.1; K4 = 0.1;
p = [tauz K1 K2 K12 K3 K4 tauw tauy];
K3 = 1; K4 = Inf;
p1 = [tauz K1 K2 K12 K3 K4 tauw tauy];


%
% Perfect adaptation FF/FB sim
%
ftnhand = @ftn_goentoro_dde;
[y0,z0,w0] = ftn_goentoro_ss(p,F0,n);
Y0 = [y0; z0; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
t = soln.x';
Y = soln.y';
z = Y(:,2);

% Then FF only sim
[y0,z0FF,w0] = ftn_goentoro_ss(p1,F0,n);
Y0 = [y0; z0FF; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p1,n,disthand);
tFF = soln.x';
Y = soln.y';
zFF = Y(:,2);


%
% Near perfect adaptation (plus), FF/FB sim
%
p(4) = K12_NPAplus;
[y0,z0plus,w0] = ftn_goentoro_ss(p,F0,n);
Y0 = [y0; z0plus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
tplus = soln.x';
Y = soln.y';
zplus = Y(:,2);

% Then FF only sim
p1(4) = K12_NPAplus;
[y0,z0FFplus,w0] = ftn_goentoro_ss(p1,F0,n);
Y0 = [y0; z0FFplus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p1,n,disthand);
tFFplus = soln.x';
Y = soln.y';
zFFplus = Y(:,2);


%
% Near perfect adaptation (plus), FF only sim
%
p(4) = K12_NPAminus;
[y0,z0minus,w0] = ftn_goentoro_ss(p,F0,n);
Y0 = [y0; z0minus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
tminus = soln.x';
Y = soln.y';
zminus = Y(:,2);

% Then FF only sim
p1(4) = K12_NPAminus;
[y0,z0FFminus,w0] = ftn_goentoro_ss(p1,F0,n);
Y0 = [y0; z0FFminus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p1,n,disthand);
tFFminus = soln.x';
Y = soln.y';
zFFminus = Y(:,2);


%
% Plotting
%
figure
plot(t,z,tplus,zplus,tminus,zminus,'linewidth',2)
set(gca,'fontsize',24)
% print(gcf,'Figs/NPA_simulation_FFFB.eps','-depsc')
% print(gcf,'Figs/NPA_simulation_FFFB.jpg','-djpeg','-r150')

figure
plot(t,z/z0,tplus,zplus/z0plus,tminus,zminus/z0minus,'linewidth',2)
hold on
plot(xlim,[1 1]*1.05,'k:',xlim,[1 1]*0.95,'k:')
set(gca,'fontsize',24)
% print(gcf,'Figs/NPA_simulation_FFFB_norm.eps','-depsc')
% print(gcf,'Figs/NPA_simulation_FFFB_norm.jpg','-djpeg','-r150')

plot(tFF,zFF/z0FF,tFFplus,zFFplus/z0FFplus,tFFminus,zFFminus/z0FFminus,'linewidth',2)
% print(gcf,'Figs/NPA_simulation_FFFB_norm1.eps','-depsc')
% print(gcf,'Figs/NPA_simulation_FFFB_norm1.jpg','-djpeg','-r150')



%}



%% ========================================================================
% Now we run a parameter sweep to get the peak value parameter P while
% varying values of K1,K2, at K12_PA, for FFL only. It turns out that this
% parameter is very weak for Region I. So that's a no-go region.
% =========================================================================
%{
nK = 50;
K11 = logspace(-2,6,nK)';
K22 = logspace(-6,2,nK)';
x0 = 1;
x1 = 10;

%
% K12_PA
%
% y0 = x0./(K11 + x0);
% y1 = x1./(K11 + x1);
% y00 = repmat(y0',nK,1);
% y11 = repmat(y1',nK,1);
% RHS = y11 - y00;
% LHS = (1 + (1./K22)*y0')/F0 - (1 + (1./K22)*y1')/F1;
% K12_PA = RHS./LHS;

%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,x0,x1);

%
% Param values
%
n = 1;
thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
% K1 = 1; K2 = 0.3; K12 = 0.3; K3 = 0.1; K4 = 0.1;
% K5 = C13*K1*K3; K6 = C23*K2*K3; K7 = C13*C23*K12*K3;
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 100];

%
% All (slow) simulation commands in the following script, which is
% commented out now, since the sims were already run and stored in the Mat
% file below.
%
% script_FigSXX_FFFB_dde_K1K2sweep_coop
load Mat/FigSXX_peak_K1K2sweep_FFFB_coop

% %
% % Boundary of Region I (and II, too)
% %
% load Mat/Fig3_FBL_K12NPA_K12sweep K1_I_FBL K2_I_FBL K1_II_FBL K2_II_FBL 



%}








