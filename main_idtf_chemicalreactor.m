%% Identifiability analysis of a chemical reactor model (DAE index 1).
% Evaluates the algebraic rank condition for identifiability over a given
% system trajectory (for a particular choice of initial conditions). 
% The DAE model of the system dynamics is described by f1 (differential
% eqs) and f2 (algebraic eqs). The measurement function is defined by h.
% The parameters sought to be identified are listed in theta.

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

%% Chemical reactor
% Step-by-step identifiability test; sensor placement.
clear all; close all; clc;

% System dimension
n1 = 2;             % differential variables
n2 = 1;             % algebraic variables
n = n1 + n2;        % full state vector

% State variables
syms t x1(t) x2(t) x3(t)
X1 = [x1(t); x2(t)];  % x1 = [c T] (reactant concentration and temperature)
X2 = [x3(t)];         % x2 = [R]   (reaction rate)
X = [X1; X2];      

% Parameters sought to be identified
syms c0 T0 Tc k1 k2 k3 k4 k5 real
                      % Tc, c0 - initial temperature, feed reactant conc.
                      % k      - constants
                      % Tc     - cooling temperature (control input)
theta = [Tc]'
p = length(theta);

% System equations
f1 = [- diff(X(1),t) + k1*(c0 - X(1)) - X(3);
      - diff(X(2),t) + k1*(T0 - X(2)) + k2*X(3) - k3*(X(2) - Tc)];
f2 =  X(3) - k5*exp(-k4/X(2))*X(1);

% Measurement function
h  = [X(1)];

% Extended system
for i = 1:p        % parameters are converted to time dependent theta(t)
    theta_t(i,1) = str2sym([sym2str(theta(i)) '(t)']);
end
X = [X1; theta_t; X2]         % augmented state vector
f1 = subs(f1,theta,theta_t);
f2 = subs(f2,theta,theta_t);
f = [f1; -diff(X(n1+1:n1+p,1)); f2]   % augmented system
h = subs(h,theta,theta_t)     % augmented measurement function
n = n + p;                    % augmented dimension

% Observability matrix
[Oc,F,H] = DAEobsvmatrix(X,f,h,n,n-2);

% Forces thetadot = thetaddot = ... = 0 in the observability matrix
xbar = [];                   % vector of symbolic state derivatives
for i = 0:size(F,2)
    xbar = [xbar; sym(['x',num2str(i)],[1 n],'real')'];
end
for i = 1:p
    Oc = subs(Oc,xbar(n1+i:n:end,1),zeros(n,1));
end

%% Chemical reactor simulation
addpath('./DAE models/')

% Chemical reactor parameters
CR.Ea = 72750; CR.k0 = 7.2e10; CR.R = 8.314; CR.V = 100; CR.rho = 1000;
CR.Cp = 0.239; CR.DeltaHr = -50e3; CR.UA = 50e3; CR.q = 100;
CR.cf = 1; CR.Tf = 350; CR.c0 = 1; CR.T0 = 350; CR.Tc = 305;

% Parameters
param.c0 = CR.c0;
param.T0 = CR.T0;
param.Tc = CR.Tc;
param.k1 = CR.q/CR.V;
param.k2 = - CR.DeltaHr/(CR.rho*CR.Cp);
param.k3 = CR.UA/(CR.V*CR.rho*CR.Cp);
param.k4 = CR.Ea/CR.R;
param.k5 = CR.k0;
param.Tc = CR.Tc;
theta_val = [param.c0 param.T0 param.Tc param.k1 param.k2 param.k3 param.k4 param.k5];

% DAE settings
x0 = [0.5; 350];        % initial value of concentration and temperature
x0(3,1) = param.k5*exp( -param.k4/x0(2) )*x0(1)      % conservation law
x0(4:n,1) = [param.Tc];

tspan = [0:0.01:10];
M = diag([1 1 0 ones(1,p)]);
options = odeset('Mass',M,'RelTol',1e-20,'AbsTol',1e-15*ones(1,length(x0)));
[t,x] = ode15s(@(t,x)chemicalreactor(t,x,theta_val,p),tspan,x0,options);
constraint1 =  param.k5*exp( -param.k4./x(:,2) ).*x(:,1);

figure(1);
subplot(131); plot(t,x(:,1)); ylabel('x1'); xlabel('t')
subplot(132); plot(t,x(:,2)); ylabel('x2'); xlabel('t')
subplot(133); plot(t,constraint1); ylabel('x3'); xlabel('t')

% Evaluate observability rank (using data)
N = length(x);     % number of data points
step = 1;          % time step between data points
T = N;             % simulation stop

% Rearrange rows in state vector to [X1 theta X2]
x = x(:,[1:n1 n-p+1:n n1+1:n1+n2]);

% Vector of symbolic derivatives
xbar = [];
for i = 0:size(F,2)
    xbar = [xbar; sym(['x',num2str(i)],[1 n],'real')'];
end

% Computes derivatives from data
xder = x;      % state
for i = 1:size(F,2)         % derxivative of order i
    for j = 1:n
        if sum(j == [n1+1:n1+p]) < 1
            xder(:,i*n+j) = gradient(xder(:,(i-1)*n+j)) ./ gradient(t(:));
        else
            xder(:,i*n+j) = zeros(N,1);
        end
    end
end

% Evaluates observability rank
count = 0;
tol = 1e-8;
for k = 1:step:T
    if mod(k,10) == 0; disp(['Counting ',num2str(k),'/',num2str(N)]); end
    count = count + 1;
    
    % Evaluates the observability rank at each time step
    Oc_k = double(subs(Oc,[xbar;c0;T0;Tc;k1;k2;k3;k4;k5],[xder(k,:)';theta_val']));

    % Identifiability
    M1_k = Oc_k(:,[n1+1:n1+p]);
    M2_k = Oc_k(:,[1:n1 n1+p+1:size(Oc,2)]);
    ident(k,1) = ( rank(Oc_k,tol) == p + rank(M2_k,tol) );
end 

%% Plot -- State evolution
figure(3)
colormap(flip(turbo))

subplot(211)
patch([t(1:step:T,1);NaN], [x(1:step:T,1);NaN], [double(ident(:)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',5,'Marker','.','MarkerSize',5,'MarkerFaceColor','flat') %
colorbar; caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_1');
box on

subplot(212)
patch([t(1:step:T,1);NaN], [x(1:step:T,2);NaN], [double(ident(:)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',5,'Marker','.','MarkerSize',5,'MarkerFaceColor','flat') %
colorbar; caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_1');
box on
