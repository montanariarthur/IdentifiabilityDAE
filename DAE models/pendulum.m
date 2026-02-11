function f = pendulum(t,x,theta,p)
% DAE model of a pendulum swing, index 3.
%
% Inputs:
%   x     -    state variables
%   theta -    parameters
%   p     -    parameters sought to be identified (thetadot = 0)
%
% References:
%
%   [1] A. N. Montanari, F. Lamoline, J. Goncalves. Identifiability of 
%       Differential-Algebraic Systems. Under review (2023).

% Copyright (C) 2023  Arthur Montanari
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% The full text of the GNU General Public License can be found in the 
% file license.txt.

% Parameters
m = theta(1);   % mass
g = theta(2);   % gravity
L = theta(3);   % length
k = theta(4);   % drag
u = theta(5);   % input

% Input

% Functions
f = [x(3);
     x(4);
     x(5)*x(1)       - k*abs(x(3))*x(3) + u;
     x(5)*x(2) - m*g - k*abs(x(4))*x(4);
     x(5)*(x(1)^2 + x(2)^2)/m + (x(3)^2 + x(4)^2) - g*x(2) + x(1)*u/m;
     zeros(p,1)];
