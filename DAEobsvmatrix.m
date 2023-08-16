function [Oc, F, H] = DAEobsvmatrix(X,f,h,n,mu,nu)
% Computes the observability matrix of a nonlinear DAE system
% F(x,xdot) = 0 with measurement function y = h(x).
%
% Inputs:
%   X     -    variables (symbolic)
%   f     -    nonlinear system (symbolic)
%   h     -    measurement function (symbolic)
%   n     -    number of state variables
%   p     -    number of parameters
%   s     -    maximum Lie derivative (default is n-1)
%
% Outputs:
%   Oc    -    observability matrix (symbolic)
%   Lieh  -    Lie derivatives (symbolic)
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

if nargin < 5      % s = n - 1 is the default input
    mu = n - 1;
    nu = n - 1;
elseif nargin < 6
    nu = mu;
end
sigma = max(mu,nu);

%% Function derivatives

% State derivatives
syms t; Xdiff = diff(X,t);

% Output derivatives
H{1} = h;
for i = 1:nu
    H{i+1} = diff(H{i},t);
end

% Function derivatives
F{1} = f;
for i = 1:mu
    F{i+1} = diff(F{i},t);
    F{i+1} = subs(F{i+1},dirac(X),zeros(n,1));  % handles dirac terms
end

%% Rename variables

% Create vector xbar = [x xdot xddot xdddot ...]
xbar = [];
for i = 0:sigma+1
    xbar = [xbar; sym(['x',num2str(i)],[1 n],'real')'];
end

% Create vector ybar = [x diff(x,t) diff(x,t,t) ...]
ybar = [X]; Xdiff = [X];
for i = 1:sigma+1
    Xdiff = diff(Xdiff, t);
    ybar = [ybar; Xdiff];
end

% Rename diff(x1, t, t, t) as x31, and so on, in functions H and F
for i = 1:nu+1
    for j = 1:length(h)
        H{i}(j) = subs(H{i}(j),flip(ybar),flip(xbar));
    end
end
for i = 1:mu+1
    for j = 1:n
        F{i}(j) = subs(F{i}(j),flip(ybar),flip(xbar));
    end
end

%% Observability matrix
M1 = []; M2 = [];
for i = 0:sigma
    M1 = vertcat(M1, F{i+1});
    M2 = vertcat(M2, H{i+1});
end
M = [M1; M2];

Oc = simplify( jacobian(M,xbar) );
% Oc = simplify( subs(Oc,dirac(xbar),zeros(length(xbar),1)) );  
% handles dirac terms (uncomment above line if needed)

end