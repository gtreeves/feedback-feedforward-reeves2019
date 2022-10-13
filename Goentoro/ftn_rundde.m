function soln = ftn_rundde(p,theta,F0,F1,n,yesplot,tspan)



if ~exist('tspan','var')
	tspan = [0 100];
end
if ~exist('yesplot','var')
	yesplot = false;
end
options = ddeset('RelTol',1e-10,'AbsTol',1e-10);

%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Initial conditions
%
[y0,z0,w0] = ftn_goentoro_ss(p,F0,n);

%
% FF/FB sim
%
ftnhand = @ftn_goentoro_dde;
Y0 = [y0; z0; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
t = soln.x';
Y = soln.y';
z = Y(:,2);

if yesplot
	plot(t,z)
end

