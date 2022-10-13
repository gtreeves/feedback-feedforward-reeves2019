function f = ftn_FB_K1K2_sweep_w_coop(p,x0,x1,n,vareps)

K1 = p(2); 
K2 = p(3);
K12 = p(4);
K3 = p(5);
C13 = p(9);
C23 = p(10);
[~,~,w0] = ftn_goentoro_ss(p,x0,n);
[~,~,w1] = ftn_goentoro_ss(p,x1,n);

y0 = x0./(K1 + x0);
y1 = x1./(K1 + x1);

N = K1*(y1*(1 + w1/C13/C23/K3) - y0/(1+vareps)*((1 + w0/C13/C23/K3)));
D = 1/(1+vareps)*(1 + K1/x0 + K1*y0/K2/x0 + w0/K3*(1/C13 + K1/x0 + K1*y0/C23/K2/x0)) - ...
	(1 + K1/x1 + K1*y1/K2/x1 + w1/K3*(1/C13 + K1/x1 + K1*y1/C23/K2/x1));

f = N/D - K12;