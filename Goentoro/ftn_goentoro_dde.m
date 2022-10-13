function dYdt = ftn_goentoro_dde(t,Y,Z,theta,p,n,disthand)
%
%
%function dYdt = ftn_goentoro_wfbk_SS(t,Y,p,F)
%
% This function contains the full equations from the Goentoro et al. model.

% % Querying the type of function calling this one:
% a = dbstack;
% if ~any(isnan(t)) && ~isempty(a) && ~any(strfind(a(2).name,'ode'))
% 	F = p;
% 	p = Y;
% 	Y = t;
% end

%
% Unpack state variables
%
y = Y(1);
z = Y(2);
w = Y(3);

%
% Unpack delayed st vbls
%
Y_thetaz = Z(:,2);
y_thetaz = Y_thetaz(1);
w_thetaz = Y_thetaz(3);
Y_thetaw = Z(:,3);
z_thetaw = Y_thetaw(2);

%
% Computing the input (delayed)
%
thetay = theta(1);
thetaz = theta(2); 
F_thetaz = disthand(t-thetaz);
F_thetay = disthand(t-thetay);


%
% Unpack parameters
%
tauz = p(1); K1 = p(2); K2 = p(3); K12 = p(4); K3 = p(5); K4 = p(6); 
tauw = p(7);
np = length(p);
if np > 7, tauy = p(8); else tauy = 1; end
if np > 8, C13 = p(9); else C13 = 1; end
if np > 9, C23 = p(10); else C23 = 1; end

%
% Param groups
%
K5 = C13*K1*K3; K6 = C23*K2*K3; K7 = C13*C23*K12*K3;

%
% Transcriptional activation functions
%
fy = (F_thetay/K1).^n./(1 + (F_thetay/K1).^n);
fz = (F_thetaz/K1).^n./(1 + (F_thetaz/K1).^n + (y_thetaz/K2).^n + ...
	(F_thetaz.*y_thetaz/K12).^n + ...
	(w_thetaz/K3).^n + (F_thetaz.*w_thetaz/K5).^n + ...
	(y_thetaz.*w_thetaz/K6).^n + ...
	(F_thetaz.*y_thetaz.*w_thetaz/K7).^n);
fw = (z_thetaw/K4).^n./(1 + (z_thetaw/K4).^n);


%
% Dynamic equations
%
dydt = (fy - y)/tauy;
dzdt = (fz - z)/tauz;
dwdt = (fw - w)/tauw;

%
% Packing the derivatives back up
%
dYdt = [dydt; dzdt; dwdt];

