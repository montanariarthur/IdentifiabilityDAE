function f = pendulum_ode(t,x,theta)
% DAE model of a pendulum swing, index 1.

% Parameters
m = theta(1);   % mass
g = theta(2);   % gravity
L = theta(3);   % length
k = theta(4);   % drag

% Functions
f = [x(3);
     x(4);
     x(5)*x(1)/m     - k*abs(x(3))*x(3)/m;
     x(5)*x(2)/m - g - k*abs(x(4))*x(4)/m;
     3*m*g/(L^2)*x(4)];
