function [y,z,w] = ftn_goentoro_ss(p,x,n)


K1 = p(2); K2 = p(3); K12 = p(4); K3 = p(5); K4 = p(6); 

np = length(p);
if np > 8
	C13 = p(9);
	if np == 9
		john
	end
else
	C13 = 1;
end
if np > 9
	C23 = p(10);
else
	C23 = 1;
end
K5 = C13*K1*K3; K6 = C23*K2*K3; K7 = C13*C23*K12*K3;

%
% Initial conditions
%
y = x.^n./(K1.^n + x.^n);

if isinf(K3) || isinf(K4)
	%
	% FFL only
	%
	z = (x/K1).^n./(1 + (x/K1).^n + (y/K2).^n + (x*y/K12).^n);
	w = 0;
	
else
	%
	% Combined FF/FB (or FB only)
	%

	if n == 1
		A0 = x./K1;
		B0 = 1 + A0 + y./K2 + x.*y./K12;
		C0 = 1/K3 + x./K5 + y./K6 + x.*y./K3./K12;
		z = -(K4.*B0-A0)/2./(B0+C0) + ...
			sqrt((K4.*B0-A0).^2 + 4*(B0+C0).*A0.*K4)/2./(B0+C0);
		w = z/K4./(1 + z/K4);
	else
		fw = @(z) (z/K4).^n./(1 + (z/K4).^n);
		fhandle = @(z) (x/K1).^n./(1 + (x/K1).^n + (y/K2).^n + ...
			(x.*y/K12).^n + (fw(z)/K3).^n + (x.*fw(z)/K5).^n + ...
			(y.*fw(z)/K6).^n + (x.*y.*fw(z)/K7).^n) - z;
		
		z = -1;
		while z < 1e-15 || z >= 1
			zguess = rand;
			z = fzero(fhandle,zguess);
		end
		w = fw(z);
	end

end