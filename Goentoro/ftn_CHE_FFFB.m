function ftn_CHE_FFFB

close all

%
% PID tuning parameters.
%
Kc = -0.1;
tauI = 1;
tauD = 0.1;

%
% System parameters
%
tspan = [-20 250];
Y01 = [50; 10]; % 50 deg C and 10 kg initially
F1 = 5; % kg/min
F2 = 5; % kg/min
% T1 is defined in our ftn
T2 = 75; % deg C

%
% Initial conditions for PID control
%
Tsp = 50; % set point (deg C)
Ts0 = 50; % initial sensor value (deg C)
c0 = 5; % F1,spec initially (kg/min)
F10 = 5; % F1 initially (kg/min)
Y0 = [Y01; Ts0; c0; F10;]; % the initial values of T and M are the same as before


%
% Other parameters
%
taua = 2/60; % 2s in minutes
taus = 6/60; % 6s in minutes


% =========================================================================
% Combined FF/FB control
% =========================================================================

%
% Initial conditions for FB/FF control (using positional form)
%
c0 = 0; % want the controller signal to start at zero b/c FF control takes care of it
I0 = 0; % integral of the error initially zero
T1s0 = 25; % initial sensor value for T1
Y0 = [Y01; Ts0; I0; F10; T1s0]; % the initial values of T and M are the same as before
ftnhand = @ftn_CHE_FFFB1;

%
% Running the integration
%
[t,Y] = ode15s(ftnhand,tspan,Y0,[],Tsp,F2,T2,Kc,tauI,tauD,taua,taus,c0);
T = Y(:,1); T = mDU(T) + 2;
T1 = 4 + 5*(t > 0);


% =========================================================================
% Plotting
% =========================================================================

%
% Plot the open and closed loop T, as well as T1
%
figure('paperpositionm','aut')
plot(t/2,T1,'Linewidth',2)
hold on
plot((t/2+2.5),T,'Linewidth',2) % added in

set(gca,'fontsize',24,'ytick',[])

xlabel('time')
xlim([-5 105])
ylim([0 10])
print(gcf,'Figs/Fig1C_CHE.jpg','-djpeg','-r150')
print(gcf,'Figs/Fig1C_CHE.eps','-depsc')









% =========================================================================
function dYdt = ftn_CHE_FFFB1(t,Y,Tsp,F2,T2,Kc,tauI,tauD,taua,taus,c0)
%
%
% This is the function for thermal mixing tank. 

%
% Unpacking our variables
%
T = Y(1);
M = Y(2);
Ts = Y(3);
I = Y(4);
F1 = Y(5);
T1s = Y(6);

%
% Evaluating our disturbance temperature at time t
%
T1 = 25 + 10*(t > 0);

%
% Defining F in terms of M
%
kprime = 3.1623;
F = kprime*sqrt(M);

%
% Sensor equations
%
dTsdt = 1/taus*(T - Ts);
dT1sdt = 1/taus*(T1 - T1s);


%
% FB Controller equation
%
e = Tsp - Ts;
c1 = c0 + Kc*(e + 1/tauI*I - tauD*dTsdt);

%
% FF controller equation
%
c2 = F2*(T2 - T)/(T - T1s);
F1spec = c1 + c2;

%
% Error integral equation
%
dIdt = e;


%
% Mass balance equation
%
dMdt = F1 + F2 - F;

%
% Energy balance equation
%
dTdt = (F1*T1 + F2*T2 - (F1 + F2)*T)/M;


%
% Actuator equation
%
dF1dt = 1/taua*(F1spec - F1);

%
% Putting them together.
%
dYdt = [dTdt; dMdt; dTsdt; dIdt; dF1dt; dT1sdt];