% script_rundde
%
% just choose your param values and the dd will run.



%
% Param values
%
n = 1;
F0 = 1;
F1 = 10;
theta = [0.5 0.5 0.5];
% K1 = 0.0212; K2 = 2.33; K3 = 0.3; K4 = Inf; K8 = Inf;
K1 = 100; K2 = 0.01; K3 = 0.3; K4 = 0.01; K8 = 0.01;
K5 = K1*K4; K6 = K2*K4; K7 = K3*K4;
tauz = 1; tauw = 1; tauy = 1;

% In case you're looking for PA
y0 = F0./(K1 + F0);
y1 = F1./(K1 + F1);
RHS = y1' - y0';
LHS = (1/F0-1/F1) - (1./K2)*(1./(K1' + F1) - 1./(K1' + F0));
K3 = RHS./LHS;


p = [tauz K1 K2 K3 K4 K8 tauw tauy];
yesplot = true;
figure
ftn_rundde(p,theta,F0,F1,n,yesplot)





% y0 = F0.^n./(K1.^n + F0.^n);
% 
% if isinf(K4) || isinf(K8)
% 	%
% 	% FFL only
% 	%
% 	z0 = (F0/K1).^n./(1 + (F0/K1).^n + (y0/K2).^n + (F0*y0/K3).^n);
% 	w0 = 0;
% 	
% else
% 	%
% 	% Combined FF/FB (or FB only)
% 	%
% 
% 	if n == 1
% 		A0 = F0./K1;
% 		B0 = 1 + A0 + y0./K2 + F0.*y0./K3;
% 		C0 = 1/K4 + F0./K5 + y0./K6 + F0.*y0./K4./K3;
% 		z0 = -(K8.*B0-A0)/2./(B0+C0) + ...
% 			sqrt((K8.*B0-A0).^2 + 4*(B0+C0).*A0.*K8)/2./(B0+C0);
% 		w0 = z0/K8./(1 + z0/K8);
% 	else
% 		fw = @(z) (z/K8).^n./(1 + (z/K8).^n);
% 		fhandle = @(z) (F0/K1).^n./(1 + (F0/K1).^n + (y0/K2).^n + ...
% 			(F0.*y0/K3).^n + (fw(z)/K4).^n + (F0.*fw(z)/K5).^n + ...
% 			(y0.*fw(z)/K6).^n + (F0.*y0.*fw(z)/K7).^n) - z;
% 		
% 		z0 = -1;
% 		while z0 < 1e-15 || z0 >= 1
% 			zguess = rand;
% 			z0 = fzero(fhandle,zguess);
% 		end
% 		w0 = fw(z0);
% 	end
% 
% end
