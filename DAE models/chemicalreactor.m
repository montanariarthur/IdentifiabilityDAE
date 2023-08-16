function f = chemicalreactor(t,x,theta,p)
% DAE model of a chemical reactor, index 1.
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
c0 = theta(1);
T0 = theta(2);
Tc = theta(3);
k(1:5) = theta(4:8);

% Functions
f = [k(1)*(c0 - x(1)) - x(3);
     k(1)*(T0 - x(2)) + k(2)*x(3) - k(3)*(x(2) - Tc);
     x(3) - k(5)*exp(-k(4)/x(2))*x(1);
     zeros(p,1)];