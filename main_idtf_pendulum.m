%% Identifiability analysis of a pendulum equation (DAE index 3).
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


%% Pendulum equation -- index 3 DAE
clear all; close all; clc;

% System dimension
n1 = 4;             % differential variables
n2 = 1;             % algebraic variables
n = n1 + n2;        % full state vector

% State variables
syms t x1(t) x2(t) x3(t) x4(t) x5(t)
X1 = [x1(t); x2(t); x3(t); x4(t)];      % x1 = [x y u v]
X2 = [x5(t)];                           % x2 = [lambda]
X = [X1; X2];      

% Parameters sought to be identified
syms m L g kdrag u real
theta = [m g L]'
p = length(theta);
kdrag = 0;

% System equations
f1 = [- diff(X(1),t) + X(3);
      - diff(X(2),t) + X(4);
      - m*diff(X(3),t) + X(5)*X(1)       - kdrag*sqrt(X(3)^2)*X(3) + u;
      - m*diff(X(4),t) + X(5)*X(2) - m*g - kdrag*sqrt(X(4)^2)*X(4)];
f2 =  X(1)^2 + X(2)^2 - L^2;

% Measurement function
h  = [atan(-X(1)/X(2))];

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
[Oc,F,H] = DAEobsvmatrix(X,f,h,n,n-1);

% Forces thetadot = thetaddot = ... = 0 in the observability matrix
xbar = [];                   % vector of symbolic state derivatives
for i = 0:size(F,2)
    xbar = [xbar; sym(['x',num2str(i)],[1 n],'real')'];
end
for i = 1:p
    Oc = subs(Oc,xbar(n1+i+n:n:end,1),zeros(size(F,2),1));
end

%% Pendulum simulation
addpath('./DAE models/')
n1 = 4;
n2 = 1;
n = n1 + n2 + p;

% Parameters
m_val = 0.3;
g_val = 9.81;
L_val = 6.25;
k_val = 0;
u_val = -1;
theta_val = [m_val g_val L_val k_val u_val];

% DAE settings
x0 = [6.25 0 0 0 -u_val/L_val m_val]';
tspan = [0:0.01:15];
M = diag([1 1 m_val m_val 0 ones(1,p)]);
options = odeset('Mass',M,'RelTol',1e-30,'AbsTol',1e-12*ones(1,length(x0)));
[t,x] = ode15s(@(t,x)pendulum(t,x,theta_val,p),tspan,x0,options);

angle = atan(-x(:,1)./x(:,2));
constraint1 = sqrt(x(:,1).^2+x(:,2).^2);
constraint2 = x(:,1).*x(:,3) + x(:,2).*x(:,4);
constraint3 = x(:,5).*(x(:,1).^2 + x(:,2).^2)/m_val + (x(:,3).^2 + x(:,4).^2) - g_val*x(:,2);

figure(1);
subplot(221); plot(t,x(:,[1 2]),t,constraint1); legend('x_1','x_2');
subplot(222); plot(t,x(:,[3 4]),t,constraint2); legend('x_3','x_4');
subplot(223); plot(t,x(:,5));              legend('x_5');
subplot(224); plot(t,angle);                    legend('y = \theta');

%% Evaluate observability rank (using data)
N = length(x);      % number of data points
step = 1;           % time step between data points
T = N;              % simulation stop

% Rearrange rows in state vector to [X1 theta X2]
x = x(:,[1:n1 n-p+1:n n1+1:n1+n2]);

% Vector of symbolic derivatives
xbar = [];
for i = 0:size(F,2)
    xbar = [xbar; sym(['x',num2str(i)],[1 n],'real')'];
end

% Computes derivatives from data
xder = x;      % state
for i = 1:size(F,2)         % derivative of order i
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
tol = 1e-6;

for k = 1:step:T
    if mod(k,10) == 0; disp(['Counting ',num2str(k),'/',num2str(N)]); end
    count = count + 1;
    
    % Evaluates the observability rank at each time step
    Oc_k = double(subs(Oc,[xbar;m;g;L;kdrag;u],[xder(k,:)';theta_val']));
    
    % Identifiability
    M1_k = Oc_k(:,[n1+1:n1+p]);
    M2_k = Oc_k(:,[1:n1 n1+p+1:size(Oc,2)]);
    ident(k,1) = ( rank(Oc_k,tol) == p + rank(M2_k,tol) )
end

%% Plot
figure(4)
colormap(flip(turbo))

subplot(121)
patch([t(1:step:T,1);NaN], [x(1:step:T,1);NaN], [double(ident(:)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',2,'Marker','o','MarkerFaceColor','flat')
colorbar; %set(gca,'ColorScale','log');
caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_2');

subplot(122)
patch([x(1:step:T,1);NaN], [x(1:step:T,2);NaN], [double(ident(:)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',2,'Marker','o','MarkerFaceColor','flat')
colorbar; %set(gca,'ColorScale','log');
caxis([-0.3 1.3]);
xlabel('x_1'); ylabel('x_2');
