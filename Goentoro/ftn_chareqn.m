function F = ftn_chareqn(X,tauz,tauw,thetaz,thetaw,J)

alpha = X(1);
beta = X(2);

falpha = tauw.*tauz.*(alpha.^2 - beta.^2) + (tauw + tauz)*alpha + ...
	1 - J.*exp(-(thetaw+thetaz).*alpha).*cos((thetaw+thetaz).*beta);
fbeta = 2*tauw.*tauz.*alpha.*beta + (tauw + tauz).*beta + ...
	J.*exp(-(thetaw+thetaz).*alpha).*sin((thetaw+thetaz).*beta);

F = [falpha; fbeta];