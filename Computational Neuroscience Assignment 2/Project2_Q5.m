clear;
clc;
close all;
fprintf('Question 5\n'); 
%% Parameters
gca = 4.4;
gk  = 8.0;
gl  = 2.0;
vca = 120;
vk  = -84;
vl  = -60;
V1 = -1.2;
V2 = 18;
V3 = 2;
V4 = 30;
C = 20;

Iext = 0;

%% Morris-Lecar functions
minf = @(V) 0.5*(1 + tanh((V-V1)/V2));
winf = @(V) 0.5*(1 + tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Equilibrium
F = @(V) -gca*minf(V).*(V-vca) ...
         -gk*winf(V).*(V-vk) ...
         -gl*(V-vl) + Iext;

Veq = fzero(F,-60);
weq = winf(Veq);

%% Simulation settings
tspan = [0 300];

% Small perturbation from equilibrium to generate the response
Y0 = [Veq + 1; weq];

%% ==============================================================
% phi = 0.01
% ==============================================================

phi = 0.01;

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;
    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

[t01,Y01] = ode45(MLE,tspan,Y0);

%% ==============================================================
% phi = 0.02
% ==============================================================

phi = 0.02;

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;
    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

[t02,Y02] = ode45(MLE,tspan,Y0);

%% ==============================================================
% phi = 0.04
% ==============================================================

phi = 0.04;

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;
    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

[t04,Y04] = ode45(MLE,tspan,Y0);

%% ==============================================================
% Action potential comparison
% ==============================================================

figure;

plot(t01,Y01(:,1),'LineWidth',1.5);
hold on;
plot(t02,Y02(:,1),'LineWidth',1.5);
plot(t04,Y04(:,1),'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');
title('Morris-Lecar Action Potential for Different \phi');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

%% ==============================================================
% Phase-plane comparison
% ==============================================================

figure;

plot(Y01(:,1),Y01(:,2),'LineWidth',1.5);
hold on;
plot(Y02(:,1),Y02(:,2),'LineWidth',1.5);
plot(Y04(:,1),Y04(:,2),'LineWidth',1.5);

plot(Veq,weq,'ko','MarkerSize',8,'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');
title('Morris-Lecar Phase Plane for Different \phi');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04','Equilibrium');
grid on;