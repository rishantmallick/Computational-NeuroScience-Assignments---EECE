clear;
clc;
close all;
fprintf('Question 6\n'); 
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
V5 = 2;
V6 = 30;
C = 20;
Iext = 0; 
%% Morris-Lecar functions
minf = @(V) 0.5*(1 + tanh((V-V1)/V2));
winf = @(V) 0.5*(1 + tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));
%% Find resting equilibrium
F = @(V) (-gca*minf(V).*(V-vca)-gk*winf(V).*(V-vk)-gl*(V-vl) + Iext);
Veq = fzero(F,-60);
weq = winf(Veq);
fprintf('Resting equilibrium:\n');
fprintf('Veq = %.6f mV\n',Veq);
fprintf('weq = %.6f\n\n',weq);
%% Morris-Lecar ODE
MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca)-gk*Y(2)*(Y(1)-vk)-gl*(Y(1)-vl) + Iext)/C;
    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];
%% ================================================================
%  INITIAL DEPOLARIZATIONS
%  ================================================================
% Initial voltages above resting potential
V0_values = Veq + [1 2 3 4 5 6 7 8 9 10 12 15 20];
% Simulation duration
tspan = [0 300];
%% Store maximum voltage
Vmax_values = zeros(size(V0_values));
%% Store trajectories
trajectories = cell(length(V0_values),1);
%% ================================================================
%  SIMULATE EACH INITIAL CONDITION
%  ================================================================
for i = 1:length(V0_values)
   V0 = V0_values(i);
   % Start w at equilibrium
   Y0 = [V0; weq];  
   % Solve MLE
   [t,Y] = ode45(MLE,tspan,Y0);
   V = Y(:,1);
   w = Y(:,2);
   % Maximum voltage reached
   Vmax_values(i) = max(V);
    % Store trajectory
    trajectories{i}.V = V;
    trajectories{i}.w = w;
    trajectories{i}.t = t;
     fprintf('V0 = %.6f mV   Vmax = %.6f mV\n', ...
            V0,Vmax_values(i));
end
figure;
hold on;

for i = 1:length(V0_values)

    plot(trajectories{i}.V,trajectories{i}.w, ...
         'LineWidth',1.2, ...
         'DisplayName',sprintf('V_0 = %.3f',V0_values(i)));
end

plot(Veq,weq,'ko','MarkerSize',8, ...
     'MarkerFaceColor','y', ...
     'LineWidth',2, ...
     'DisplayName','Equilibrium');

xlabel('V (mV)');
ylabel('w');
title('Morris-Lecar Phase-Plane Trajectories');
legend('Location','best');
grid on;

%% ================================================================
%  MAXIMUM AMPLITUDE VS INITIAL VOLTAGE
%  ================================================================

figure;

plot(V0_values,Vmax_values,'ko-', ...
     'LineWidth',1.5, ...
     'MarkerSize',5);

xlabel('Initial voltage V_0 (mV)');
ylabel('Maximum voltage V_{max} (mV)');
title('Maximum Action-Potential Amplitude vs Initial Voltage');
grid on;
V0_values = Veq + (0:60);
% Simulation duration
tspan = [0 300];
%% Store maximum voltage
Vthreshold_values = zeros(size(V0_values));
%% Store trajectories
trajectories = cell(length(V0_values),1);
%% ================================================================
%  SIMULATE EACH INITIAL CONDITION
%  ================================================================
for i = 1:length(V0_values)
   V0 = V0_values(i);
   % Start w at equilibrium
   Y0 = [V0; weq];  
   % Solve MLE
   [t,Y] = ode45(MLE,tspan,Y0);
   V = Y(:,1);
   w = Y(:,2);
   % Maximum voltage reached
   Vthreshold_values(i) = max(V);
    % Store trajectory
    trajectories{i}.V = V;
    trajectories{i}.w = w;
    trajectories{i}.t = t;
     fprintf('V0 = %.6f mV   Vmax = %.6f mV\n', ...
            V0,Vthreshold_values(i));
end
figure;

plot(V0_values,Vthreshold_values,'ko-', ...
     'LineWidth',1.5, ...
     'MarkerSize',5);

xlabel('Initial voltage V_0 (mV)');
ylabel('Maximum threshold voltage V_{threshold max} (mV)');
title('Maximum Action threshold -Potential Amplitude vs Initial Voltage');
grid on;
V0_below = Veq + 40;
V0_above = Veq + 50;

figure;
hold on;

% Below threshold
[t1,Y1] = ode45(MLE,tspan,[V0_below;weq]);

% Above threshold
[t2,Y2] = ode45(MLE,tspan,[V0_above;weq]);

plot(Y1(:,1),Y1(:,2),'LineWidth',2, ...
     'DisplayName',sprintf('Below: V_0 = %.3f',V0_below));

plot(Y2(:,1),Y2(:,2),'LineWidth',2, ...
     'DisplayName',sprintf('Above: V_0 = %.3f',V0_above));

plot(Veq,weq,'ko','MarkerSize',8, ...
     'MarkerFaceColor','y', ...
     'DisplayName','Equilibrium');

xlabel('V (mV)');
ylabel('w');
title('Phase Plane: Below and Above Threshold');
legend('Location','best');
grid on;