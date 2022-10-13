% function f = ftn_RegionI_K1upper_w_coop(p,x0,x1,n,vareps)
function f = ftn_RegionI_K1upper_w_coop(x,K3,K4,C13,C23,x0,x1,n,vareps)

K1 = exp(x);
p = [1 K1 Inf Inf K3 K4 1 1 C13 C23];
% K1 = p(2);
% K3 = p(5);
% C13 = p(8);

[~,~,w0] = ftn_goentoro_ss(p,x0,n);
[~,~,w1] = ftn_goentoro_ss(p,x1,n);

N = 1/(1+vareps)*(1 + w0/C13/K3) - (1 + w1/C13/K3);
D = 1/x1*(1 + w1/K3) - 1/(1+vareps)/x0*(1 + w0/K3);

f = N/D - K1;