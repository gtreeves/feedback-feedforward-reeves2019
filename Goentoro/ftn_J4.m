function [g,J4,J5] = ftn_J4(Ki,K1,K2,K3,K4,Kname,J,F,n,FF,Y0,K12)
%
%
% This function evaluates J4, which is the partial of fz wrt w, evaluated
% at steady state. All four parameters (K1,K2,K3,K4) are needed to evaluate
% the steady state. F is the fold-change from basal F0 = 1, and n is the
% hill coefficient.
% 
% The logical switch "FF" will tell us whether FF control is present.
%
% There are some strange input parameters, so let me explain. I want to
% keep this function as modular as possible, so in theory we could fix any
% three of the four input "K's" and vary the third (using fzero) to find
% where g = 0 = J4*J5 - J.
%
% "Ki": the varied parameter to try to get g = 0, if using fzero
% "K1" thru "K4": obvious
% "Kname": the name of the varied parameter. So if we are holding K1,K2,K4
%	fixed, and varying K3, then "Kname" would be = 'K3' (a character
%	variable).
% "J": the value of J4*J5 that we are targeting. This comes from an
%	analysis of the characteristic equation.
% "F": the fold-change F
% "n": Hill coeff
% "FF": logical true if we are looking at the feedforward model.

eval([Kname,' = ',num2str(Ki),';'])
if ~exist('K12','var')
	K12 = K1*K2; 
end

%
% Steady state
%
fw = @(z) (z/K4).^n./(1 + (z/K4).^n);
if ~exist('Y0','var')
	p = [1 K1 K2 K12 K3 K4];
	[y0,z0,w0] = ftn_goentoro_ss(p,F,n);
else
	y0 = Y0(1);
	z0 = Y0(2);
	w0 = Y0(3);
end

%
% Next, define fz as a function of w (fw as a function of z was already
% defined)
%
K5 = K1*K3; K6 = K2*K3; K7 = K12*K3;

% if FF % FF/FB 
% 	fz = @(w) (F/K1).^n./(1 + (F/K1).^n + (y0/K2).^n + ...
% 			(F.*y0/K12).^n + (w/K3).^n + (F.*w/K5).^n + ...
% 			(y0.*w/K6).^n + (F.*y0.*w/K7).^n);
% else % FB only
% 	fz = @(w)(F/K1).^n./(1 + (F/K1).^n + (w/K3).^n + (F.*w/K5).^n);
% end
A = (F/K1).^n;
if FF % FF/FB
	B = 1 + (F/K1).^n + (y0/K2).^n + ...
			(F.*y0/K12).^n;
	C = (1/K3).^n + (F/K5).^n + ...
			(y0/K6).^n + (F.*y0/K7).^n;
else % FB only
	B = 1 + (F/K1).^n;
	C = (1/K3).^n + (F/K5).^n;
end
fz = @(w)A/(B + C*w.^n);



%
% Now, calculate J4
%
% delt = 1e-4;
% fz1 = fz(w0);
% fz2 = fz(w0*(1 + delt));
% J41 = (fz2 - fz1)/w0/delt;

J4 = -A*C*n*w0.^(n-1)./(B + C*w0.^n).^2;

%
% Finally, J5
%
% fw1 = fw(z0);
% fw2 = fw(z0*(1 + delt));
% J51 = (fw2 - fw1)/z0/delt;

a = z0/K4;
J5 = (n*a.^(n-1))./(1 + a.^n).^2/K4;

%
% Now, pack them up:
%
g = J4.*J5 - J;













