% script_Topic08_mixing_tank_soln
%
% This script is for a FF control example from Topic09

close all

%
% PID tuning parameters.
%
Kc = -1;
tauI = 1;
tauD = 0.1;

%
% System parameters
%
tspan = [0 30];
Y01 = [50; 10]; % 50 deg C and 10 kg initially
F1 = 5; % kg/min
F2 = 5; % kg/min
% T1 is defined in our ftn
T2 = 75; % deg C
ftnhand = @ftn_Topic08_mixing_tank_soln;

%
% Integrate the open loop using ode15s.
%
[tOL,Y1] = ode15s(ftnhand,tspan,Y01,[],F1,F2,T2);

% unpack our output as one of the elements of ode15s's output
TOL = Y1(:,1);
MOL = Y1(:,2);

% =========================================================================
% Now doing PID control
% =========================================================================

%
% Initial conditions for PID control
%
Tsp = 50; % set point (deg C)
Ts0 = 50; % initial sensor value (deg C)
c0 = 5; % F1,spec initially (kg/min)
F10 = 5; % F1 initially (kg/min)
Y0 = [Y01; Ts0; c0; F10;]; % the initial values of T and M are the same as before
ftnhand = @ftn_Topic08_mixing_tank_PID_soln;

%
% Other parameters
%
taua = 2/60; % 2s in minutes
taus = 6/60; % 6s in minutes

%
% Running the integration
%
[t,Y] = ode15s(ftnhand,tspan,Y0,[],Tsp,F2,T2,Kc,tauI,tauD,taua,taus);

%
% Unpacking state variables
%
T = Y(:,1);
M = Y(:,2); % controller signal
Ts = Y(:,3);
c = Y(:,4);
F1 = Y(:,5);

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
ftnhand = @ftn_Topic09_mixing_tank_FFFB_soln;

%
% Running the integration
%
[tFF,Y] = ode15s(ftnhand,tspan,Y0,[],Tsp,F2,T2,Kc,tauI,tauD,taua,taus,c0);

%
% Unpacking state variables
%
TFF = Y(:,1);
MFF = Y(:,2); 
TsFF = Y(:,3);
IFF = Y(:,4); % integral of the error
F1FF = Y(:,5);
T1sFF = Y(:,6);

%
% Controller signal
%
e = Tsp - TsFF;
dTsdt = 1/taus*(TFF - TsFF);
c1 = c0 + Kc*(e + 1/tauI*IFF - tauD*dTsdt);
c2 = F2*(T2 - TFF)./(TFF - T1sFF);
F1spec = c1 + c2;


% =========================================================================
% Plotting
% =========================================================================

%
% First we plot the open and closed loop T, as well as T1
%
P = 10; % minutes
T1 = 25 + 10*sin(pi*t/P).*(t < P); % defining T1
figure
plot(t,T,'b','Linewidth',2)
hold on
plot(tOL,TOL,'g','Linewidth',2)
plot(tFF,TFF,'c','Linewidth',2) % added in
plot(t,T1,'r -.','Linewidth',2)

xlabel('time [min]')
ylabel('Temperature [^oC]')
legend('Closed loop outlet T','Open loop outlet T','FF/FB T','Disturbance T1') % changed
ylim([20 60])

%
% Next we plot just the closed loop T to see it better
%
figure
plot(t,T,'Linewidth',2) % what is the temperature really doing?
hold on
plot(tFF,TFF,'c','Linewidth',2) % added in

xlabel('time [min]')
ylabel('Temperature [^oC]')
legend('Closed loop outlet T','FF/FB T') % changed

%
% Finally we plot the controller signals
%
figure
plot(t,c,'b','Linewidth',2)
hold on
plot(tFF,c1,'g','Linewidth',2)
plot(tFF,c2,'c','Linewidth',2)
plot(tFF,F1spec,'r -.','Linewidth',2)

xlabel('time [min]')
ylabel('Controller signals')
legend('Closed loop c','FF/FB c1','FF/FB c2','FF/FB F1spec')



% %
% % Finally we plot flow rate F1 and M.
% %
% figure
% plot(t,M,'b','Linewidth',2)
% hold on
% plot(tOL,MOL,'g','Linewidth',2)
% plot(tFF,MFF,'c','Linewidth',2)
% plot(t,F1,'r -.','Linewidth',2) % trade-offs: you gain mass
% 
% xlabel('time [min]')
% ylabel('Tank mass [kg] or inlet flow [kg/min]')
% legend('Closed loop mass','Open loop mass','FF/FB Mass','Inlet flow F1')







