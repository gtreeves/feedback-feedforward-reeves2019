function D = ftn_unitstep(t,t0,x0,x1)
%
%
%
%
% This function provides a description of the disturbance variable.  Here
% we provide it as a unit step.

D = zeros(size(t));
D(t < t0) = x0;
D(t >= t0) = x1;