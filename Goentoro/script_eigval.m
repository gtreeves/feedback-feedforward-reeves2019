% script_eigval

clear
close all

F1 = 10;
K1 = 1; 
K2 = 0.3;
K4 = 0.1;
K8 = 0.1;
n = 2;

%
% In this first section, we acquire the beta that gives us the instability
% point. This "betastar" is the frequency of oscillations at instability.
% Note that this beta (and subsequently, J) is ONLY dependent on tauw,
% tauz, thetaw, and thetaz. That is an interesting observation, because the
% steady states are independent of those parameters.
%
tauw = 1;
tauz = 1;
thetaz = 0.5;
thetaw = 0.5;

betamax = pi/2/(thetaw + thetaz) - 0.05;
beta = linspace(0,betamax)';
fbeta = eigval_alpha0(beta,tauz,tauw,thetaz,thetaw);

betastar = advNtnDU(@eigval_alpha0,betamax/2,[1e-4,betamax],[],[],[],[],[],...
	tauz,tauw,thetaz,thetaw);

%
% Now that we have beta, we calculate J, which is the other varible
% parameter in the characteristic equation. In a very complicated way, all
% other parameters (K1...K8, as well as the steady state values of our
% variables) are embedded in J.
%
J = (1-tauw*tauz*betastar.^2)./cos((thetaw+thetaz).*betastar);



%% ========================================================================
% Now that we have J, we can find the stability boundaries for different
% scenarios. In this scenario, we will vary K2 and find the boundary for
% K4.
% =========================================================================
%{


%
% Given J, can we reverse engineer the relationships between the other
% model parameters? J = J4*J5, where J4 is the partial of fz wrt w, and J5
% is the partial of fw wrt z, (both evaluated at ss).
%
% Here we run fzero on a function that gives us J4*J5 - J as an output,
% where J4 and J5 are the partials evaluated inside the function, and J is
% the input paramter, whose value is given by the above equation. We will
% st K1 and K8 to given values (at top of script), then sweep K2 and find
% the K4 that gives us J4*J5 - J = 0 
%
nK2 = 300;
K22 = flipud(logspace(-2,0,nK2)');
K44 = zeros(nK2,1);

Kname = 'K4';
FF = true;
fhandle = @ftn_J4;
Ki = 0.08431;
K4_0 = 1; % placeholder

count = 1;
while Ki > 1e-2
	K44(count) = fzero(fhandle,Ki,[],K1,K22(count),K4_0,K8,Kname,J,F1,n,FF);
	Ki = K44(count);
	count = count + 1;
end
if count <= nK2
	K22(count:end) = [];
	K44(count:end) = [];
end
[K22,isort] = sort(K22);
K44 = K44(isort);

%
% Now we find the K4 that gives instability for our chosen F1, K1, K8, n
% parameters in the FB only scenario.
%
FF = false;
K2_0 = 1; % placeholder...means nothing in FB only
K44_FB = fzero(fhandle,Ki,[],K1,K2_0,K4_0,K8,Kname,J,F1,n,FF);

figure
loglog(K22,K44)
hold on
plot([1e-2,1],K44_FB*[1 1])
xlim([1e-2 1])
ylim([1e-2 1])
xlabel('K_2')
ylabel('K_4')
legend('FF/FB','FB only')

%}




%% ========================================================================
% In this scenario, we will vary K8 and find the boundary for K4.
% =========================================================================
%{


nK8 = 300;
K88 = logspace(-2,0,nK8)';
K44 = zeros(nK8,1);

Kname = 'K4';
FF = true;
fhandle = @ftn_J4;
Ki = 0.08;
K4_0 = 1; % placeholder

count = 1;
while Ki > 1e-2 && count <= nK8
	K44(count) = fzero(fhandle,Ki,[],K1,K2,K4_0,K88(count),Kname,J,F1,n,FF);
	Ki = K44(count);
	count = count + 1;
end
if count <= nK8
	K88(count:end) = [];
	K44(count:end) = [];
end

%
% Now we find the K4's that give instability for our chosen F1, K1, K2, n
% parameters in the FB only scenario.
%
K88_FB = logspace(-2,0,nK8)';
K44_FB = zeros(nK8,1);

FF = false;
Ki = 0.03;
K2_0 = 1; % placeholder...means nothing in FB only
K4_0 = 1; % placeholder

count = 1;
while Ki > 1e-2 && count <= nK8
	K44_FB(count) = fzero(fhandle,Ki,[],K1,K2,K4_0,K88_FB(count),Kname,J,F1,n,FF);
	Ki = K44_FB(count);
	count = count + 1;
end
if count <= nK8
	K88_FB(count:end) = [];
	K44_FB(count:end) = [];
end


figure
loglog(K88,K44)
hold on
loglog(K88_FB,K44_FB)
xlim([1e-2 1])
ylim([1e-2 1])
xlabel('K_8')
ylabel('K_4')
legend('FF/FB','FB only')

%}




%% ========================================================================
% Here we will vary beta from betastar down to zero and trace out the
% "nullclines" in the complex plane.
% =========================================================================
% {

nb = 500;
Beta = linspace(betastar,0,nb)';
dbeta = (Beta(2) - Beta(1));
fhandle = @ftn_nullcline_real;
alpha0 = -1e-4;
Alpha1 = zeros(nb,1);
for i = 1:nb

% 	Alpha(i) = advNtnDU(fhandle,alpha0,[-100 0],[],[],[],[],[],...
% 		Beta(i),tauz,tauw,thetaz,thetaw,J);
	Alpha1(i) = fzero(fhandle,alpha0,[],Beta(i),tauz,tauw,thetaz,thetaw,J);
	if isnan(Alpha1(i))
		Beta(i:end) = [];
		Alpha1(i:end) = [];
		break
	end
	if i == 1
		alpha0 = Alpha1(i);
	else
		m = (Alpha1(i) - Alpha1(i-1))./(Beta(i) - Beta(i-1));
		alpha0 = alpha0 + m*dbeta;
	end
end


fhandle = @ftn_nullcline_imag;
alpha0 = -1e-4;
Alpha2 = zeros(nb,1);
for i = 1:nb

% 	Alpha(i) = advNtnDU(fhandle,alpha0,[-100 0],[],[],[],[],[],...
% 		Beta(i),tauz,tauw,thetaz,thetaw,J);
	Alpha2(i) = fzero(fhandle,alpha0,[],Beta(i),tauz,tauw,thetaz,thetaw,J);
	if isnan(Alpha2(i))
		Beta(i:end) = [];
		Alpha2(i:end) = [];
		break
	end
	if i == 1
		alpha0 = Alpha2(i);
	else
		m = (Alpha2(i) - Alpha2(i-1))./(Beta(i) - Beta(i-1));
		alpha0 = alpha0 + m*dbeta;
	end
end


figure
plot(Alpha1,Beta)
hold on
plot(Alpha2,Beta)




%}


%%
% falpha = fhandle(linspace(-10,0)',Beta(i),tauz,tauw,thetaz,thetaw,J);
% a = tauw*tauz;
% b = tauw + tauz;
% c = 1 - tauw*tauz.*Beta.^2;
% p = @(alpha,beta)a.*alpha.^2 + b.*alpha + 1 - tauw*tauz.*beta.^2;
% g = @(alpha,beta)J.*exp(-(thetaw+thetaz).*alpha).*cos((thetaw+thetaz).*beta);
% 
% Alpha = linspace(-2,0)';
% figure
% plot(Alpha,[p(Alpha,Beta(30)) g(Alpha,Beta(30))])

nb = 50;
Beta = linspace(betastar,0,nb)';
p = @(alpha,beta)2*tauw.*tauz.*alpha.*beta + (tauw + tauz).*beta;
g = @(alpha,beta)-J.*exp(-(thetaw+thetaz).*alpha).*sin((thetaw+thetaz).*beta);

Alpha = linspace(-2,0)';
figure
plot(Alpha,[p(Alpha,Beta(2)) g(Alpha,Beta(2))])


% %
% % OK, that's not working...fzero is never finding a root. I hope I didn't
% % do some algebra wrong, because there's no check point before this.
% % Anyway, we'll brute force it by calculating J4, J5 for many K4's at K2 =
% % 0.2.
% %
% K2 = 0.5;
% Kname = 'K4';
% FF = true;
% fhandle = @ftn_J4;
% Ki = 0.5;
% K4_0 = 1; % placeholder
% K44 = logspace(-2,0)';
% for i = 1:length(K44)
% 	[g(i),J4(i),J5(i)] = ftn_J4(K44(i),K1,K2,K4_0,K8,Kname,J,F1,n,FF);
% end













