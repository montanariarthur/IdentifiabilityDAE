%% Identifiability analysis of a linear DAE model.
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


%% Identifiability of linear DAE system - order 4
% Which structures are identifiable?
clear all; close all; clc;

% System dimension
n1 = 2;             % differential variables
n2 = 2;             % algebraic variables
n = n1 + n2;        % full state vector

% State variables
syms x1(t) x2(t) x3(t) x4(t)
X1 = [x1(t);x2(t)];
X2 = [x3(t);x4(t)];
X = [X1; X2];

% Parameters sought to be identified
A = sym('A',[n n])
A11 = A(1:n1,1:n1);       A12 = A(1:n1,n1+1:n1+n2);
A21 = A(n1+1:n1+n2,1:n1); A22 = A(n1+1:n1+n2,n1+1:n1+n2);

theta = A(:);
p = length(theta);

% System equations
E = blkdiag(eye(n1),zeros(n2));
f = - E*diff(X,t) + A*X;

% Measurement function
C = diag([1 1 1 1]);
h = C*X

% Extended system
for i = 1:p        % parameters are converted to time dependent theta(t)
    theta_t(i,1) = str2sym([sym2str(theta(i)) '(t)']);
end
X = [X1; theta_t; X2]         % augmented state vector
f = subs(f,theta,theta_t);
f = [f(1:n1,1); -diff(X(n1+1:n1+p,1)); f(n1+1:n1+n2,1)]   % augmented
h = subs(h,theta,theta_t)     % augmented measurement function
n = n + p;                    % augmented dimension

% Observability matrix
[Oc,F,H] = DAEobsvmatrix(X,f,h,n,6-1);

% Forces thetadot = thetaddot = ... = 0 in the observability matrix
xbar = [];                   % vector of symbolic state derivatives
for i = 0:size(F,2)
    xbar = [xbar; sym(['x',num2str(i)],[1 n],'real')'];
end
for i = 1:p
    Oc = subs(Oc,xbar(n1+i+n:n:end,1),zeros(size(F,2),1));
end

% Rank condition
M1 = Oc(:,[n1+1:n1+p]);
M2 = Oc(:,[1:n1 n1+p+1:size(Oc,2)]);
rankOc = rank(Oc);
rankM2 = rank(M2);
if rankOc == p + rank(M2) 
    disp(['Augmented DAE system is identifiable.'])
else
    disp(['Augmented DAE system is *not* identifiable.'])
end

%% Linear DAE simulation
addpath('./DAE models/')
n1 = 2;
n2 = 2;
n = n1 + n2 + p;

% Parameters
sys = rss(n1+n2);
A_val = sys.A;
A11_val = A_val(1:n1,1:n1);       A12_val = A_val(1:n1,n1+1:n1+n2);
A21_val = A_val(n1+1:n1+n2,1:n1); A22_val = A_val(n1+1:n1+n2,n1+1:n1+n2);

% DAE settings
x01 = randn(n1,1);
x02 = -inv(A22_val)*A21_val*x01;
theta_val = A_val(:); theta_val = theta_val(theta_val~=0);
x0 = [x01; x02; theta_val];

tspan = [0:0.01:10];
M = blkdiag(E,eye(p));
options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',1e-10*ones(1,length(x0)));
[t,x] = ode15s(@(t,x)linearDAE(t,x,A_val,p),tspan,x0,options);

figure(1);
subplot(121); plot(t,x(:,[1 2])); ylabel('x1');
subplot(122); plot(t,x(:,[3 4])); ylabel('x2');


%% Evaluate observability rank (using data)
N = length(x);      % number of data points
step = 10;          % time step between data points
tol = 1e-6;
T = N;           % simulation stop

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
tol_range = 10.^(-[1:15]);

for k = 1:step:T
    if mod(k,1) == 0; disp(['Counting ',num2str(k),'/',num2str(N)]); end
    count = count + 1;
    
    % Evaluates the observability rank at each time step
    Oc_k = double(subs(Oc,[xbar],[xder(k,:)']));
    
    % Identifiability
    M1_k = Oc_k(:,[n1+1:n1+p]);
    M2_k = Oc_k(:,[1:n1 n1+p+1:size(Oc,2)]);
    ident(k,1) = ( rank(Oc_k,tol) == p + rank(M2_k,tol) );
end

%% Plot
figure(4)
colormap(flip(turbo))

subplot(221)
patch([t(1:step:T,1);NaN], [x(1:step:T,1);NaN], [double(ident(1:step:end)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',2,'Marker','o','MarkerFaceColor','flat')
colorbar; %set(gca,'ColorScale','log');
caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_2');

subplot(222)
patch([t(1:step:T,1);NaN], [x(1:step:T,2);NaN], [double(ident(1:step:end)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',2,'Marker','o','MarkerFaceColor','flat')
colorbar; %set(gca,'ColorScale','log');
caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_2');

subplot(223)
patch([t(1:step:T,1);NaN], [x(1:step:T,3);NaN], [double(ident(1:step:end)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',2,'Marker','o','MarkerFaceColor','flat')
colorbar; %set(gca,'ColorScale','log');
caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_3');

subplot(224)
patch([t(1:step:T,1);NaN], [x(1:step:T,4);NaN], [double(ident(1:step:end)); double(ident(end))],...
    'EdgeColor','flat','LineWidth',2,'Marker','o','MarkerFaceColor','flat')
colorbar; %set(gca,'ColorScale','log');
caxis([-0.3 1.3]);
xlabel('t'); ylabel('x_4');
