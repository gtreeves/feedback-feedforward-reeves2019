function f = ftn_RegionI_II_FB_w_coop(x,K2,K3,K4,C13,C23,x0,x1,n,vareps)


K1 = exp(x);
p = [1 K1 K2 Inf K3 K4 1 1 C13 C23];
[~,~,w0] = ftn_goentoro_ss(p,x0,n);
[~,~,w1] = ftn_goentoro_ss(p,x1,n);

y0 = x0./(K1 + x0);
y1 = x1./(K1 + x1);

% A = (1 + w0/C23/K3)./(1 + w1/C23/K3);
% N = A/(1+vareps)*(K1*y0/x0) - (K1*y1/x1);

% D = 1/(1+vareps)/(1 + w1/C23/K3)*((1 + K1/x0 + w0/K3*(1/C13 + K1/x0)) - ...
% 	(1 + K1/x1 + w1/K3*(1/C13 + K1/x1)));

% D = (1 + K1/x1 + w1/K3*(1/C13 + K1/x1))/(1 + w1/C23/K3) - ...
% 	1/(1+vareps)/(1 + w1/C23/K3)*(1 + K1/x0 + w0/K3*(1/C13 + K1/x0));


N = (1 + w0/C23/K3)/(1+vareps)*K1*y0/x0 - (1 + w1/C23/K3)*K1*y1/x1;
D = (1 + K1/x1 + w1/K3*(1/C13 + K1/x1)) - ...
	1/(1+vareps)*(1 + K1/x0 + w0/K3*(1/C13 + K1/x0));

f = N/D - K2;