function fbeta = eigval_alpha0(beta,tauz,tauw,thetaz,thetaw)
%
% This function is the imaginary part of the characteristic equation after
% the real part is solved for J and plugged in

fbeta = (tauw + tauz)*beta + ...
	(1-tauw*tauz*beta.^2).*tan((thetaw+thetaz).*beta);