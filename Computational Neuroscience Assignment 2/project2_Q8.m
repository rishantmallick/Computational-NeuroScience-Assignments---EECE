clear;
clc;
close all;
fprintf('Question 8\n'); 
%% Parameters
gca = 4.4;
gk  = 8.0;
gl  = 2.0;
vca = 120;
vk  = -84;
vl  = -60;
phi = 0.02;

V1 = -1.2;
V2 = 18;
V3 = 2;
V4 = 30;
C = 20;

Iext = 86;

%% Morris-Lecar functions
minf = @(V) 0.5*(1 + tanh((V-V1)/V2));
winf = @(V) 0.5*(1 + tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% MLE equations
MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;
     
    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

%% ---------------------------------------------------------------
% Find equilibrium for Iext = 86
% ---------------------------------------------------------------

F = @(V) -gca*minf(V).*(V-vca) ...
         -gk*winf(V).*(V-vk) ...
         -gl*(V-vl) + Iext;

Veq = fzero(F,-30);
weq = winf(Veq);

fprintf('Equilibrium:\n');
fprintf('Veq = %.6f mV\n',Veq);
fprintf('weq = %.6f\n',weq);

%% ---------------------------------------------------------------
% Run the model BACKWARDS in time
% ---------------------------------------------------------------

% Choose an initial point away from the equilibrium
Y0 = [-27.9; 0.17];

% Negative time means backwards integration
tspan_backward = [0 -5000];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[tb,Yb] = ode45(MLE,tspan_backward,Y0,options);

%% ---------------------------------------------------------------
% Plot backward trajectory
% ---------------------------------------------------------------

figure;

plot(Yb(:,1),Yb(:,2),'LineWidth',1.2);
hold on;

plot(Veq,weq,'ko','MarkerSize',8,...
     'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');
title('Backward-Time Trajectory');
legend('Backward trajectory','Equilibrium');
grid on;

%% ---------------------------------------------------------------
% Zoom into the periodic orbit
% ---------------------------------------------------------------

figure;

plot(Yb(:,1),Yb(:,2),'LineWidth',1.5);
hold on;

plot(Veq,weq,'ko','MarkerSize',8,...
     'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');
title('Unstable Periodic Orbit (UPO)');
legend('UPO','Equilibrium');
grid on;