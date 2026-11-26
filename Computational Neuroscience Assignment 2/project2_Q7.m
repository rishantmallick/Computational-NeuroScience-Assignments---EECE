clear;
clc;
close all;
fprintf('Question 7\n'); 
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

%% ================================================================
% Equilibrium for Iext = 0
% Used for trial 1
% ================================================================

I0 = 0;

F0 = @(V) -gca*minf(V).*(V-vca) ...
          -gk*winf(V).*(V-vk) ...
          -gl*(V-vl) + I0;

Veq0 = fzero(F0,-60);
weq0 = winf(Veq0);

fprintf('Equilibrium for Iext = 0:\n');
fprintf('Veq0 = %.6f mV\n',Veq0);
fprintf('weq0 = %.6f\n\n',weq0);

%% ================================================================
% Equilibrium for Iext = 86
% Used for trial 2
% ================================================================

F86 = @(V) -gca*minf(V).*(V-vca) ...
           -gk*winf(V).*(V-vk) ...
           -gl*(V-vl) + Iext;

Veq86 = fzero(F86,-30);
weq86 = winf(Veq86);

fprintf('Equilibrium for Iext = 86:\n');
fprintf('Veq86 = %.6f mV\n',Veq86);
fprintf('weq86 = %.6f\n\n',weq86);

%% ================================================================
% Morris-Lecar equations
% ================================================================

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;
     
    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

%% Simulation time
tspan = [0 1000];

%% ================================================================
% THREE INITIAL CONDITIONS
% ================================================================

% Trial 1: equilibrium for Iext = 0
Y01 = [Veq0; weq0];

% Trial 2: equilibrium for Iext = 86
Y02 = [Veq86; weq86];

% Trial 3: specified off-equilibrium point
Y03 = [-27.9; 0.17];

%% ================================================================
% SIMULATE
% ================================================================

[t1,Y1] = ode45(MLE,tspan,Y01);
[t2,Y2] = ode45(MLE,tspan,Y02);
[t3,Y3] = ode45(MLE,tspan,Y03);

%% ================================================================
% PHASE PLANE
% ================================================================

figure;
hold on;

plot(Y1(:,1),Y1(:,2),'LineWidth',1.5);
plot(Y2(:,1),Y2(:,2),'LineWidth',1.5);
plot(Y3(:,1),Y3(:,2),'LineWidth',1.5);

% Equilibrium for Iext = 86
plot(Veq86,weq86,'ko','MarkerSize',8, ...
     'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');

title('Morris-Lecar Phase Plane, I_{ext} = 86 \muA/cm^2');

legend( ...
    'Trial 1: I_{ext}=0 equilibrium', ...
    'Trial 2: I_{ext}=86 equilibrium', ...
    'Trial 3: (-27.9, 0.17)', ...
    'Equilibrium for I_{ext}=86', ...
    'Location','best');

grid on;