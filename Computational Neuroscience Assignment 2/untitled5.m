clear;
clc;
close all;
fprintf('Question 2 and 3\n'); 
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

% Anonymous Functions
minf = @(V) 0.5 * (1 + tanh((V - V1) / V2));
winf = @(V) 0.5 * (1 + tanh((V - V3) / V4));
% Added + eps to prevent division by zero exactly at vk (-84 mV)
V_null = @(V) (-gca * minf(V) .* (V - vca) - gl * (V - vl) + Iext) ./ (gk * (V - vk) + eps);
tauw = @(V) 1.0 ./ cosh((V - V3) / (2 * V4));

%% ========================================================================
%  EQUILIBRIUM SOLVERS(Using NullClines)
%  =======================================================================
F = @(V) (-gca * minf(V) .* (V - vca) - gk * winf(V) .* (V - vk) - gl * (V - vl) + Iext);
Veq = fzero(F, -60);
weq = winf(Veq);

fprintf('Morris-Lecar Equilibrium Point:\n');
fprintf('Veq = %.6f mV\n', Veq);
fprintf('weq = %.6f\n\n', weq);
%% ========================================================================
%  Using fsolve to solve
%  =======================================================================
% Simultaneous equilibrium equations
f = @(x) [(-gca*minf(x(1))*(x(1)-vca)-gk*x(2)*(x(1)-vk)-gl*(x(1)-vl) + Iext)/C;
        phi*(winf(x(1))-x(2))/tauw(x(1))
];

options = optimoptions('fsolve','Display','off');

eq = fsolve(f, [-60; 0.01], options);

Veq = eq(1);
weq = eq(2);
fprintf('Solved with simultaneous equilibrium:\n');
fprintf('Veq = %.6f mV\n', Veq);
fprintf('weq = %.6f\n', weq);
%% ========================================================================
%  GRID & VECTOR CALCULATIONS
%  =======================================================================
% Dense arrays for continuous nullcline plotting curves
V_curves = linspace(-100, 100, 500);
W_w = winf(V_curves);
W_V = V_null(V_curves);

% Narrower, focused grid for the direction arrows so they don't crowd
% Spacing out steps (e.g., step size of 8 for V) lets arrow heads breath
[Vgrid, Wgrid] = meshgrid(-90:8:50, -0.1:0.08:1.1);

m = minf(Vgrid);
w_inf = winf(Vgrid);
tau = tauw(Vgrid);

dV = (-gca * m .* (Vgrid - vca) - gk * Wgrid .* (Vgrid - vk) - gl * (Vgrid - vl) + Iext) / C;
dw = phi * (w_inf - Wgrid) ./ tau;

% Pure, raw normalization without any artificial multiplier distortion
L = sqrt(dV.^2 + dw.^2);
L(L == 0) = 1;
dV_plot = dV ./ L;
dw_plot = dw ./ L;

%% ========================================================================
%  FIGURE 1: Pure Nullclines and Equilibrium Only
%  =======================================================================
figure('Name', 'Figure 1: Nullclines Only');
plot(V_curves, W_V, 'r-', 'LineWidth', 2.5, 'DisplayName', 'V-nullcline');
hold on;
plot(V_curves, W_w, 'b-', 'LineWidth', 2.5, 'DisplayName', 'w-nullcline');
plot(Veq, weq, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'y', 'LineWidth', 2, 'DisplayName', 'Equilibrium');

xlabel('V (mV)'); ylabel('w');
title('Morris-Lecar Nullclines');
legend('Location', 'best'); grid on;
xlim([-90 50]); ylim([-0.1 1.1]);

%% ========================================================================
%  FIGURE 2: Direction Field Only
%  =======================================================================
figure('Name', 'Figure 2: Direction Field Only');
% Autoscale parameter (0.4) controls arrow length relative to grid spacing
quiver(Vgrid, Wgrid, dV_plot, dw_plot, 0.4, 'k', 'LineWidth', 1.0);

xlabel('V (mV)'); ylabel('w');
title('Morris-Lecar Vector Flow Field');
grid on;
xlim([-90 50]); ylim([-0.1 1.1]);

%% ========================================================================
%  FIGURE 3: Clean, Combined Phase Plane
%  =======================================================================
figure('Name', 'Figure 3: Clean Combined Phase Plane');

% 1. Plot spacing-optimized direction field arrows in background
quiver(Vgrid, Wgrid, dV_plot, dw_plot, 0.4, 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
hold on;

% 2. Plot Nullcline Curves
plot(V_curves, W_V, 'r-', 'LineWidth', 2.5, 'DisplayName', 'V-nullcline');
plot(V_curves, W_w, 'b-', 'LineWidth', 2.5, 'DisplayName', 'w-nullcline');

% 3. Plot Equilibrium Node Point
plot(Veq, weq, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'LineWidth', 2, ...
     'DisplayName', sprintf('Equilibrium: (%.2f, %.2f)', Veq, weq));

xlabel('V (mV)'); ylabel('w');
title('Morris-Lecar Phase Plane Context');
legend('Location', 'best'); grid on;

% Force plot limits to focus right where the dynamics intersect
xlim([-90 50]); 
ylim([-0.1 1.1]);
%% ========================================================================
%  JACOBIAN AND EIGENVALUES AT EQUILIBRIUM
%  ========================================================================

% Define symbolic variables
syms V w

% Symbolic versions of the gating functions
m_sym = 0.5 * (1 + tanh((V - V1) / V2));
winf_sym = 0.5 * (1 + tanh((V - V3) / V4));
tauw_sym = 1.0 / cosh((V - V3) / (2 * V4));

% Morris-Lecar equations
f1 = (-gca*m_sym*(V-vca)-gk*w*(V-vk)-gl*(V-vl) + Iext) / C;

f2 = phi*(winf_sym-w) / tauw_sym;

% Obtain Jacobian
J = jacobian([f1; f2], [V w]);

% Evaluate Jacobian at equilibrium
J_eq = double(subs(J, [V w], [Veq weq]));

% Calculate eigenvalues
eigenvalues = eig(J_eq);

% Display results
disp('Jacobian at equilibrium:')
format long
disp(J_eq)

disp('Eigenvalues:')
disp(eigenvalues)

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


fprintf('Question 9\n'); 
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
I = 80:100; 
firingrate = zeros(size(I)); 
minf = @(V) 0.5*(1+tanh((V-V1)/V2)); 
winf = @(V) 0.5*(1+tanh((V-V3)/V4)); 
tauw = @(V) 1./cosh((V-V3)/(2*V4)); 

for i = 1:length(I) 
  %% MLE equations 
  Iext = I(i); 
  MLE = @(t,Y) [ ... 
    (-gca*minf(Y(1))*(Y(1)-vca) ... 
     -gk*Y(2)*(Y(1)-vk) ... 
     -gl*(Y(1)-vl) + Iext)/C; 
      
    phi*(winf(Y(1))-Y(2))/tauw(Y(1)) 
  ]; 

  F = @(V) -gca*minf(V).*(V-vca) ... 
         -gk*winf(V).*(V-vk) ... 
         -gl*(V-vl) + Iext; 
 
  Veq = fzero(F,-30); 
  weq = winf(Veq);  

  % Morris-Lecar equations 
  syms V w 
  f1 = (-gca*minf(V)*(V-vca)-gk*w*(V-vk)-gl*(V-vl) + Iext) / C; 
  f2 = phi*(winf(V)-w) / tauw(V); 

  % Obtain Jacobian 
  J = jacobian([f1; f2], [V w]); 

  % Evaluate Jacobian at equilibrium 
  J_eq = double(subs(J, [V w], [Veq weq])); 

  % Calculate eigenvalues 
  eigenvalues = eig(J_eq); 

  % Display results 
  fprintf('Jacobian at equilibrium: %g\n', Iext) 
  format long 
  disp(J_eq) 
  fprintf('Eigen Values at equilibrium: %g\n', Iext) 
  disp(eigenvalues) 

  tspan = [0 3000]; 
  Y0=[Veq,weq]; 
  options = odeset('RelTol',1e-8,'AbsTol',1e-10); 
  [t,Y] = ode45(MLE,tspan,Y0,options); 
  V1 = Y(:,1); 

  % 2. Find spike times 
  figure; 
  plot(t,V1); 
  xlabel('Time (ms)'); 
  ylabel('Voltage (mV)'); 
 
  titleName = ['V vs t, Iext = ', num2str(Iext)]; 
  title(titleName); 
  grid on; 
 
  filename = ['V_vs_t_Iext_' num2str(Iext) '_Q9.png']; 
  exportgraphics(gcf, filename, 'Resolution', 300); 

  [pks,locs] = findpeaks(V1,t,'MinPeakHeight',-50); 

  % 3. Calculate inter-spike intervals (ISI) 
  ISI = diff(locs); 

  % 4. Calculate firing rate 
  if ~isempty(ISI) 
    firing_rate = 1000 / mean(ISI); 
  else 
    firing_rate = 0; 
  end 

  firingrate(i) = firing_rate; 

  % Display results 
  disp('Spike times (ms):'); 
  disp(locs); 
  disp('ISI (ms):'); 
  disp(ISI); 
  fprintf('Firing rate = %.4f Hz\n',firing_rate); 
end 

figure; 
plot(I, firingrate, 'o-'); 
xlabel('I_{ext} (\muA/cm^2)'); 
ylabel('Firing Rate (Hz)'); 
title('Firing Rate vs Applied Current'); 
grid on; 

exportgraphics(gcf, 'Firing_Rate_vs_Applied_Current_Q9.png', 'Resolution', 300);


fprintf('Question 10\n');
%% Parameters
gca = 4;
gk  = 8.0;
gl  = 2;
vca = 120;
vk  = -84;
vl  = -60;
phi = 0.0667;
V1 = -1.2;
V2 = 18;
V3 = 12;
V4 = 17.4;
V5 = 12;
V6 = 17.4;
C = 20;
Iext = 30;
%% Morris-Lecar functions
minf = @(V) 0.5*(1 + tanh((V-V1)/V2));
winf = @(V) 0.5*(1 + tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));
%% Equilibrium condition
% At equilibrium:
% dV/dt = 0
% dw/dt = 0
%
% Therefore w = winf(V)
F = @(V) (-gca*minf(V).*(V-vca) ...
          -gk*winf(V).*(V-vk) ...
          -gl*(V-vl) + Iext);
%% Find the three equilibrium voltages
Vscan = linspace(-100,100,2000);
Fscan = F(Vscan);
rootsV = [];
for i = 1:length(Vscan)-1
    if Fscan(i)*Fscan(i+1)<0
        root = fzero(F,[Vscan(i),Vscan(i+1)]);
        if isempty(rootsV) || all(abs(rootsV-root) > 1e-5)
            rootsV(end+1) = root;
        end
    end
end
Veq = rootsV;
weq = winf(Veq);
disp('Equilibrium points:')
disp([Veq(:), weq(:)])
for k = 1:length(Veq)
   syms Vsym wsym
   msym = 0.5*(1 + tanh((Vsym-V1)/V2));
   winfsym = 0.5*(1 + tanh((Vsym-V3)/V4));
   tauwsym = 1/cosh((Vsym-V3)/(2*V4));
   f1 = (-gca*msym*(Vsym-vca) ...
      -gk*wsym*(Vsym-vk) ...
      -gl*(Vsym-vl) + Iext)/C;
   f2 = phi*(winfsym-wsym)/tauwsym;
   J = jacobian([f1; f2],[Vsym wsym]);
   J_eq = double(subs(J,{Vsym,wsym},{Veq(k),weq(k)}));
   lambda = eig(J_eq);
   fprintf('Equilibrium %d:\n',k);
   fprintf('Veq = %.4f, weq = %.4f\n',Veq(k),weq(k));
   disp('Jacobian:');
   disp(J_eq);
   disp('Eigenvalues:');
   disp(lambda);
   if all(real(lambda) < 0)
    disp('Stable')
   elseif all(real(lambda) > 0)
    disp('Unstable')
   elseif any(real(lambda) > 0) && any(real(lambda) < 0)
    disp('Saddle point')
   end
end

V = linspace(-100,100,2000);
% V-nullcline: dV/dt = 0
V_null = (-gca*minf(V).*(V-vca) ...
          -gl*(V-vl) + Iext) ./ (gk*(V-vk));

% w-nullcline: dw/dt = 0
w_null = winf(V);
for k = 1:length(Veq)

    J_eq = double(subs(J,{Vsym,wsym},{Veq(k),weq(k)}));
    lambda = eig(J_eq);

    if real(lambda(1))*real(lambda(2)) < 0
        saddle_index = k;
    end
end

Vs = Veq(saddle_index);
Ws = weq(saddle_index);
J_saddle = double(subs(J,{Vsym,wsym},{Vs,Ws}));
[VEC,D] = eig(J_saddle);
lambda = diag(D);
unstable_index = find(real(lambda) > 0);
stable_index   = find(real(lambda) < 0);

v_unstable = VEC(:,unstable_index);
v_stable   = VEC(:,stable_index);
odefun = @(t,y) [
    (-gca*minf(y(1))*(y(1)-vca) ...
     -gk*y(2)*(y(1)-vk) ...
     -gl*(y(1)-vl) + Iext)/C;

    phi*(winf(y(1))-y(2))/tauw(y(1))
];
epsilon = 1e-3;
y0_u1 = [Vs; Ws] + epsilon*v_unstable;
y0_u2 = [Vs; Ws] - epsilon*v_unstable;

[t_u1,Y_u1] = ode45(odefun,[0 50],y0_u1);
[t_u2,Y_u2] = ode45(odefun,[0 50],y0_u2);
y0_s1 = [Vs; Ws] + epsilon*v_stable;
y0_s2 = [Vs; Ws] - epsilon*v_stable;

[t_s1,Y_s1] = ode45(odefun,[0 -50],y0_s1);
[t_s2,Y_s2] = ode45(odefun,[0 -50],y0_s2);
%% Figure 1: Nullclines and equilibrium points

figure;
hold on;

plot(V,V_null,'LineWidth',2);
plot(V,w_null,'LineWidth',2);

plot(Veq,weq,'ko','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('w');
title('Morris-Lecar Nullclines and Equilibrium Points');

legend('V-nullcline', ...
       'w-nullcline', ...
       'Equilibrium points');

grid on;
xlim([-50 15]);
ylim([0 0.4]);

hold off;


%% Figure 2: Stable and unstable manifolds

figure;
hold on;

% Unstable manifold
hU = plot(Y_u1(:,1),Y_u1(:,2),'-','LineWidth',2.5);
plot(Y_u2(:,1),Y_u2(:,2),'-','LineWidth',2.5);

% Stable manifold
hS = plot(Y_s1(:,1),Y_s1(:,2),'--','LineWidth',2.5);
plot(Y_s2(:,1),Y_s2(:,2),'--','LineWidth',2.5);

% Saddle point
hE = plot(Vs,Ws,'ko','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('w');
title('Stable and Unstable Manifolds');

legend([hU hS hE], ...
       'Unstable manifold', ...
       'Stable manifold', ...
       'Saddle point');

grid on;
xlim([-50 15]);
ylim([0 0.4]);

hold off;

fprintf('Question 11\n');
%% Parameters
gca = 4;
gk  = 8.0;
gl  = 2;
vca = 120;
vk  = -84;
vl  = -60;
phi = 0.0667;
V1 = -1.2;
V2 = 18;
V3 = 12;
V4 = 17.4;
V5 = 12;
V6 = 17.4;
C = 20;
Iext = 30;
%% Morris-Lecar functions
minf = @(V) 0.5*(1 + tanh((V-V1)/V2));
winf = @(V) 0.5*(1 + tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));
Vscan = linspace(-100,100,2000);
I_values = linspace(30,50,200);
Veq_all = NaN(length(I_values),3);
weq_all = NaN(length(I_values),3);
for j = 1:length(I_values)
    Iext = I_values(j);
    F = @(V) (-gca*minf(V).*(V-vca) ...
          -gk*winf(V).*(V-vk) ...
          -gl*(V-vl) + Iext);
    Fscan = F(Vscan);
    rootsV = [];
    for i = 1:length(Vscan)-1
        if Fscan(i)*Fscan(i+1) < 0
          root = fzero(F,[Vscan(i),Vscan(i+1)]);
          if isempty(rootsV) || all(abs(rootsV-root)>1e-5)
                rootsV(end+1) = root;
          end
        end
    end
    Veq_all(j,1:length(rootsV)) = rootsV;
    weq_all(j,1:length(rootsV)) = winf(rootsV);
end
Vscan1 = linspace(-100,100,2000);
I_values1 = linspace(39,40,100);
Veq_all1 = NaN(length(I_values1),3);
weq_all1 = NaN(length(I_values1),3);
for j = 1:length(I_values1)
    Iext = I_values1(j);
    F = @(V) (-gca*minf(V).*(V-vca) ...
          -gk*winf(V).*(V-vk) ...
          -gl*(V-vl) + Iext);
    Fscan = F(Vscan1);
    rootsV1 = [];
    for i = 1:length(Vscan1)-1
        if Fscan(i)*Fscan(i+1) < 0
          root = fzero(F,[Vscan1(i),Vscan1(i+1)]);
          if isempty(rootsV1) || all(abs(rootsV1-root)>1e-5)
                rootsV1(end+1) = root;
          end
        end
    end
    Veq_all1(j,1:length(rootsV1)) = rootsV1;
    weq_all1(j,1:length(rootsV1)) = winf(rootsV1);
end
figure;
hold on;

plot(I_values1,Veq_all1(:,1),'LineWidth',2);
plot(I_values1,Veq_all1(:,2),'LineWidth',2);
plot(I_values1,Veq_all1(:,3),'LineWidth',2);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('V_{eq} (mV)');

title('Equilibrium Points for I_{ext} = 39 to 40 \muA/cm^2');

legend('Equilibrium 1', ...
       'Equilibrium 2', ...
       'Equilibrium 3');

grid on;

xlim([39 40]);

hold off;
%% Morris-Lecar functions

minf = @(V) 0.5*(1 + tanh((V-V1)/V2));

winf = @(V) 0.5*(1 + tanh((V-V3)/V4));

tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Current range

I_values = 30:0.25:45;

firing_rate = zeros(size(I_values));

%% Simulation settings

tspan = [0 2000];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%% Loop over current

for j = 1:length(I_values)

    Iext = I_values(j);

    %% Morris-Lecar equations

    MLE = @(t,y) [
        (-gca*minf(y(1))*(y(1)-vca) ...
        -gk*y(2)*(y(1)-vk) ...
        -gl*(y(1)-vl) + Iext)/C;

        phi*(winf(y(1))-y(2))/tauw(y(1))
    ];

    %% Initial condition

    % Start from the low-voltage equilibrium region
    Y0 = [-41.8452, 0.002047];

    %% Solve

    [t,Y] = ode45(MLE,tspan,Y0,options);

    Vtrace = Y(:,1);

    %% Find spikes

    [pks,locs] = findpeaks(Vtrace,t,...
                           'MinPeakHeight',-50);

    %% Inter-spike intervals

    ISI = diff(locs);

    %% Firing rate

    if ~isempty(ISI)

        firing_rate(j) = 1000/mean(ISI);

    else

        firing_rate(j) = 0;

    end

end

%% Display results

disp('Iext       Firing Rate (Hz)');

disp([I_values(:), firing_rate(:)]);

%% Plot firing rate vs applied current

figure;

plot(I_values,firing_rate,'o-','LineWidth',2);

xlabel('Applied current I_{ext} (\muA/cm^2)');
ylabel('Firing rate (Hz)');

title('Firing Rate vs Applied Current');

grid on;


fprintf('Question 12\n');
%% Classical Hodgkin-Huxley Model

%% Parameters

C = 1;                  % uF/cm^2

GK  = 36;               % mS/cm^2
GNa = 120;              % mS/cm^2
GL  = 0.3;              % mS/cm^2

EK  = -72;              % mV
ENa = 55;               % mV

%% Temperature

T = 6.3;
Q = 3;

phi = Q^((T-6.3)/10);

% At T = 6.3:
% phi = 1

%% External current

Iext = 10;              % uA/cm^2

%% Leak reversal potential

EL = -49.4;             % mV

%% Rate functions

% alpha_n
alpha_n = @(V) ...
    -0.01*phi*(V+50) ./ ...
    (exp(-(V+50)/10)-1);

% beta_n
beta_n = @(V) ...
    0.125*phi*exp(-(V+60)/80);

% alpha_m
alpha_m = @(V) ...
    -0.1*phi*(V+35) ./ ...
    (exp(-(V+35)/10)-1);

% beta_m
beta_m = @(V) ...
    4*phi*exp(-(V+60)/18);

% alpha_h
alpha_h = @(V) ...
    0.07*phi*exp(-(V+60)/20);

% beta_h
beta_h = @(V) ...
    phi./(exp(-(V+30)/10)+1);


%% Initial voltage

V0 = -60;


%% Handle 0/0 cases

% alpha_n has 0/0 at V = -50 mV
%
% lim V->-50 alpha_n(V) = phi

if abs(V0 + 50) < 1e-10
    an0 = 0.1* phi;
else
    an0 = alpha_n(V0);
end
% alpha_m has 0/0 at V = -35 mV
%
% lim V->-35 alpha_m(V) = phi

if abs(V0 + 35) < 1e-10
    am0 = phi;
else
    am0 = alpha_m(V0);
end
ah0 = alpha_h(V0);

bn0 = beta_n(V0);
bm0 = beta_m(V0);
bh0 = beta_h(V0);
n0 = an0/(an0 + bn0);
m0 = am0/(am0 + bm0);
h0 = ah0/(ah0 + bh0);
%% Initial condition

% Y(1) = V
% Y(2) = n
% Y(3) = m
% Y(4) = h

Y0 = [V0; n0; m0; h0];
%% HH differential equations

HH = @(t,Y) [
        % dV/dt

    (Iext ...
    - GK*Y(2)^4*(Y(1)-EK) ...
    - GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    - GL*(Y(1)-EL))/C;


    % dn/dt

    phi*( ...
    alpha_n_safe(Y(1),phi)*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2) );


    % dm/dt

    phi*( ...
    alpha_m_safe(Y(1),phi)*(1-Y(3)) ...
    - beta_m(Y(1))*Y(3) );


    % dh/dt

    phi*( ...
    alpha_h(Y(1))*(1-Y(4)) ...
    - beta_h(Y(1))*Y(4) )

];


%% Simulation time

tspan = [0 1000];


%% ODE solver options

options = odeset('RelTol',1e-8,'AbsTol',1e-10);


%% Solve the system

[t,Y] = ode45(HH,tspan,Y0,options);


%% Extract variables

V = Y(:,1);
n = Y(:,2);
m = Y(:,3);
h = Y(:,4);


%% Calculate conductances

GK_t = GK*n.^4;

GNa_t = GNa*m.^3.*h;

GL_t = GL*ones(size(t));


%% Plot membrane potential

figure;

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');

title(['Hodgkin-Huxley Model, I_{ext} = ', ...
       num2str(Iext), ' \muA/cm^2']);

grid on;


%% Plot gating variables

figure;

plot(t,n,'LineWidth',1.5);
hold on;

plot(t,m,'LineWidth',1.5);
plot(t,h,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Gating variables');

title('HH Gating Variables');

legend('n','m','h');

grid on;

hold off;


%% Plot conductances

figure;

plot(t,GK_t,'LineWidth',1.5);
hold on;

plot(t,GNa_t,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Conductance (mS/cm^2)');

title('HH Conductances');

legend('G_K','G_{Na}');

grid on;

hold off;


%% Display values

fprintf('Temperature = %.1f degrees C\n',T);
fprintf('Phi = %.4f\n',phi);

fprintf('\nInitial conditions:\n');
fprintf('V0 = %.4f mV\n',V0);
fprintf('n0 = %.4f\n',n0);
fprintf('m0 = %.4f\n',m0);
fprintf('h0 = %.4f\n',h0);
%% Safe alpha_n

function a = alpha_n_safe(V,phi)

    if abs(V + 50) < 1e-8

        % Limit at V = -50 mV
        a = 0.1 * phi;

    else

        a = -0.01*phi*(V+50) / ...
            (exp(-(V+50)/10)-1);

    end

end


%% Safe alpha_m

function a = alpha_m_safe(V,phi)

    if abs(V + 35) < 1e-8

        % Limit at V = -35 mV
        a = phi;

    else

        a = -0.1*phi*(V+35) / ...
            (exp(-(V+35)/10)-1);

    end

end



fprintf('Question 13\n');
%% Parameters

C = 1;                  % uF/cm^2

GK  = 36;               % mS/cm^2
GNa = 120;              % mS/cm^2
GL  = 0.3;              % mS/cm^2

EK  = -72;              % mV
ENa = 55;               % mV

%% Temperature

T = 6.3;
Q = 3;

phi = Q^((T-6.3)/10);

%% Determine EL for Vrest = -60 mV

Vrest = -60;
Irest = 0;

% Steady-state gating variables

ninf = alpha_n_safe(Vrest,phi) / ...
       (alpha_n_safe(Vrest,phi) + beta_n(Vrest));

minf = alpha_m_safe(Vrest,phi) / ...
       (alpha_m_safe(Vrest,phi) + beta_m(Vrest));

hinf = alpha_h(Vrest) / ...
       (alpha_h(Vrest) + beta_h(Vrest));

% Calculate EL

EL = Vrest - ...
     (-GK*ninf^4*(Vrest-EK) ...
      -GNa*minf^3*hinf*(Vrest-ENa) ...
      + Irest)/GL;

fprintf('Temperature = %.2f C\n',T);
fprintf('phi = %.4f\n',phi);

fprintf('\nResting-state gating variables:\n');
fprintf('n_inf = %.6f\n',ninf);
fprintf('m_inf = %.6f\n',minf);
fprintf('h_inf = %.6f\n',hinf);

fprintf('\nCalculated leak reversal potential:\n');
fprintf('EL = %.6f mV\n',EL);


%% Now apply external current

Iext = 10;

%% Initial condition

% Start exactly at the resting state

Y0 = [Vrest; ninf; minf; hinf];


%% Hodgkin-Huxley equations

HH = @(t,Y) [

    % dV/dt

    (Iext ...
    - GK*Y(2)^4*(Y(1)-EK) ...
    - GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    - GL*(Y(1)-EL))/C;


    % dn/dt

    phi*( ...
    alpha_n_safe(Y(1),phi)*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2));


    % dm/dt

    phi*( ...
    alpha_m_safe(Y(1),phi)*(1-Y(3)) ...
    - beta_m(Y(1))*Y(3));


    % dh/dt

    phi*( ...
    alpha_h(Y(1))*(1-Y(4)) ...
    - beta_h(Y(1))*Y(4))
];


%% Simulation

tspan = [0 1000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);


%% Use stiff solver ode15s

[t,Y] = ode15s(HH,tspan,Y0,options);


%% Extract variables

V = Y(:,1);
n = Y(:,2);
m = Y(:,3);
h = Y(:,4);


%% Plot membrane potential

figure;

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');

title(['Hodgkin-Huxley Model, I_{ext} = ', ...
       num2str(Iext), ' \muA/cm^2']);

grid on;


%% Plot gating variables

figure;

plot(t,n,'LineWidth',1.5);
hold on;

plot(t,m,'LineWidth',1.5);
plot(t,h,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Gating variable');

title('HH Gating Variables');

legend('n','m','h');

grid on;

hold off;


%% Safe alpha_n function

function a = alpha_n_safe(V,phi)

    if abs(V + 50) < 1e-8

        % Limit as V -> -50
        a = 0.1*phi;

    else

        a = -0.01*phi*(V+50) / ...
            (exp(-(V+50)/10)-1);

    end

end


%% Safe alpha_m function

function a = alpha_m_safe(V,phi)

    if abs(V + 35) < 1e-8

        % Limit as V -> -35
        a = phi;

    else

        a = -0.1*phi*(V+35) / ...
            (exp(-(V+35)/10)-1);

    end

end


%% beta_n

function b = beta_n(V)

    b = 0.125*exp(-(V+60)/80);

end


%% beta_m

function b = beta_m(V)

    b = 4*exp(-(V+60)/18);

end


%% alpha_h

function a = alpha_h(V)

    a = 0.07*exp(-(V+60)/20);

end


%% beta_h

function b = beta_h(V)

    b = 1./(exp(-(V+30)/10)+1);

end


fprintf('Question 14\n');
%% Parameters

C = 1;                  % uF/cm^2

GK  = 36;               % mS/cm^2
GNa = 120;              % mS/cm^2
GL  = 0.3;              % mS/cm^2

EK  = -72;              % mV
ENa = 55;               % mV
T = 6.3;
Q = 3;
EL = -49.4;
phi = Q^((T-6.3)/10);

%% Rate functions

alpha_n = @(V) alpha_n_safe(V,phi);

beta_n = @(V) ...
    0.125*phi.*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V,phi);

beta_m = @(V) ...
    4*phi.*exp(-(V+60)/18);

alpha_h = @(V) ...
    0.07*phi.*exp(-(V+60)/20);

beta_h = @(V) ...
    phi./(exp(-(V+30)/10)+1);

%% Steady-state gating variables

n_inf = @(V) alpha_n(V)./(alpha_n(V)+beta_n(V));

m_inf = @(V) alpha_m(V)./(alpha_m(V)+beta_m(V));

h_inf = @(V) alpha_h(V)./(alpha_h(V)+beta_h(V));


%% =========================================================
% PART 1: RESTING EQUILIBRIUM
% ==========================================================

Iext = 0;

% At equilibrium:
%
% n = n_inf(V)
% m = m_inf(V)
% h = h_inf(V)

F = @(V) ...
    -GK*n_inf(V).^4.*(V-EK) ...
    -GNa*m_inf(V).^3.*h_inf(V).*(V-ENa) ...
    -GL*(V-EL) + Iext;


%% Find equilibrium

Vscan = linspace(-100,100,3000);

Fscan = F(Vscan);

rootsV = [];

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1) < 0

        root = fzero(F,[Vscan(i),Vscan(i+1)]);

        if isempty(rootsV) || ...
                all(abs(rootsV-root) > 1e-5)

            rootsV(end+1) = root;

        end

    end

end


Vrest = rootsV(1);

nrest = n_inf(Vrest);
mrest = m_inf(Vrest);
hrest = h_inf(Vrest);


fprintf('\n====================================\n');
fprintf('RESTING EQUILIBRIUM\n');
fprintf('====================================\n');

fprintf('Vrest = %.6f mV\n',Vrest);
fprintf('nrest = %.6f\n',nrest);
fprintf('mrest = %.6f\n',mrest);
fprintf('hrest = %.6f\n',hrest);


%% =========================================================
% PART 2: JACOBIAN AND STABILITY
% ==========================================================
syms Vs ns ms hs

f1 = (-GK*ns^4*(Vs-EK) ...
      -GNa*ms^3*hs*(Vs-ENa) ...
      -GL*(Vs-EL) + Iext)/C;

alpha_n_sym = -0.01*phi*(Vs+50) / ...
              (exp(-(Vs+50)/10)-1);

beta_n_sym = 0.125*phi*exp(-(Vs+60)/80);

alpha_m_sym = -0.1*phi*(Vs+35) / ...
              (exp(-(Vs+35)/10)-1);

beta_m_sym = 4*phi*exp(-(Vs+60)/18);

f2 = phi*(alpha_n_sym*(1-ns) ...
      - beta_n_sym*ns);

f3 = phi*(alpha_m_sym*(1-ms) ...
      - beta_m_sym*ms);

f4 = phi*(alpha_h(Vs)*(1-hs) ...
      - beta_h(Vs)*hs);


J = jacobian([f1;f2;f3;f4], ...
             [Vs ns ms hs]);


Jrest = double(subs(J, ...
    {Vs,ns,ms,hs}, ...
    {Vrest,nrest,mrest,hrest}));


eigenvalues = eig(Jrest);


fprintf('\n====================================\n');
fprintf('STABILITY AT REST\n');
fprintf('====================================\n');

disp('Jacobian:');
disp(Jrest);

disp('Eigenvalues:');
disp(eigenvalues);


if all(real(eigenvalues) < 0)

    fprintf('Resting equilibrium is STABLE.\n');

elseif any(real(eigenvalues) > 0)

    fprintf('Resting equilibrium is UNSTABLE.\n');

else

    fprintf('Resting equilibrium is NON-HYPERBOLIC.\n');

end


%% =========================================================
% PART 3: HH DIFFERENTIAL EQUATIONS
% ==========================================================

HH = @(t,Y) [

    (-GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    -GL*(Y(1)-EL) + Iext)/C;

    phi*(alpha_n(Y(1))*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2));

    phi*(alpha_m(Y(1))*(1-Y(3)) ...
    - beta_m(Y(1))*Y(3));

    phi*(alpha_h(Y(1))*(1-Y(4)) ...
    - beta_h(Y(1))*Y(4))
];


%% Numerical tolerances

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);
%% =========================================================
% PART 4: THRESHOLD DEPOLARIZATION
% ==========================================================

% Brief current pulse is approximated by directly
% depolarizing the membrane.
deltaV_values = 0:0.1:40;
threshold_found = false;
spike_result = zeros(size(deltaV_values));
for k = 1:length(deltaV_values)

    deltaV = deltaV_values(k);

    % Initial condition:
    %
    % V is depolarized
    % n, m and h remain at their resting values

    Y0 = [Vrest + deltaV;
          nrest;
          mrest;
          hrest];


    [t,Y] = ode15s(HH,[0 100],Y0,options);

    V = Y(:,1);
        % Detect action potential

    if max(V) > 0

        spike_result(k) = 1;
         if ~threshold_found

            threshold_found = true;

            threshold_deltaV = deltaV;

            threshold_voltage = Vrest + deltaV;

            t_threshold = t;

            V_threshold = V;

        end

    end

end
if threshold_found

    fprintf('Threshold depolarization = %.2f mV\n', ...
        threshold_deltaV);

    fprintf('Threshold voltage = %.2f mV\n', ...
        threshold_voltage);

else

    fprintf('No spike found in the tested range.\n');

end
%% =========================================================
% Plot threshold search
% ==========================================================

figure;

plot(deltaV_values,spike_result,'o-','LineWidth',1.5);

xlabel('\DeltaV (mV)');
ylabel('Spike');

title('Threshold Depolarization');

grid on;
%% =========================================================
% Plot voltage at threshold
% ==========================================================

if threshold_found

    figure;

    plot(t_threshold,V_threshold,'LineWidth',1.5);

    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');

    title(['Response at threshold, \DeltaV = ', ...
        num2str(threshold_deltaV),' mV']);

    grid on;

end
%% =========================================================
% Safe alpha_n
% ==========================================================

function a = alpha_n_safe(V,phi)

    a = -0.01*phi.*(V+50) ./ ...
        (exp(-(V+50)/10)-1);

    idx = abs(V+50) < 1e-8;

    a(idx) = 0.1*phi;

end


%% =========================================================
% Safe alpha_m
% ==========================================================

function a = alpha_m_safe(V,phi)

    a = -0.1*phi.*(V+35) ./ ...
        (exp(-(V+35)/10)-1);

    idx = abs(V+35) < 1e-8;

    a(idx) = phi;

end



fprintf('Question 15\n');
%% Parameters

gca = 4;
gk  = 8.0;
gl  = 2;

vca = 120;
vk  = -84;
vl  = -60;

phi = 0.0667;

V1 = -1.2;
V2 = 18;

V3 = 12;
V4 = 17.4;

V5 = 12;
V6 = 17.4;

C = 20;


%% Morris-Lecar functions

minf = @(V) 0.5*(1 + tanh((V-V1)/V2));

winf = @(V) 0.5*(1 + tanh((V-V3)/V4));

tauw = @(V) 1./cosh((V-V3)/(2*V4));


%% Current range

I_values = linspace(8,12,200);

% Voltage range used to search for equilibria
Vscan = linspace(-100,100,3000);

% Maximum number of equilibrium points
max_eq = 3;

Veq_all = NaN(length(I_values),max_eq);
weq_all = NaN(length(I_values),max_eq);

% Stability:
%  1  = stable
% -1  = saddle
%  0  = unstable/non-hyperbolic

stability_all = NaN(length(I_values),max_eq);


%% Symbolic Jacobian

syms V w

m = 0.5*(1 + tanh((V-V1)/V2));

w_inf = 0.5*(1 + tanh((V-V3)/V4));

tau = 1/cosh((V-V3)/(2*V4));

%% Morris-Lecar equations

f1 = (-gca*m*(V-vca) ...
      -gk*w*(V-vk) ...
      -gl*(V-vl))/C;

f2 = phi*(w_inf-w)/tau;

J = jacobian([f1;f2],[V w]);


%% Find equilibrium points for every current

for j = 1:length(I_values)

    Iext = I_values(j);

    % Equilibrium condition:
    % dV/dt = 0 and w = winf(V)

    F = @(V) ...
        -gca*minf(V).*(V-vca) ...
        -gk*winf(V).*(V-vk) ...
        -gl*(V-vl) + Iext;

    Fscan = F(Vscan);

    rootsV = [];

    %% Search for all roots

    for i = 1:length(Vscan)-1

        if Fscan(i)*Fscan(i+1) < 0

            root = fzero(F,[Vscan(i),Vscan(i+1)]);

            if isempty(rootsV) || ...
                    all(abs(rootsV-root) > 1e-5)

                rootsV(end+1) = root;

            end
        end
    end


    %% Calculate stability of every equilibrium

    for r = 1:length(rootsV)

        Veq = rootsV(r);
        weq = winf(Veq);

        % Evaluate Jacobian at equilibrium
        J_eq = double(subs(J,{V,w},{Veq,weq}));

        % Eigenvalues
        lambda = eig(J_eq);

        % Store equilibrium
        Veq_all(j,r) = Veq;
        weq_all(j,r) = weq;


        % Stability classification

        if all(real(lambda) < 0)

            stability_all(j,r) = 1;       % Stable

        elseif any(real(lambda) > 0) && ...
               any(real(lambda) < 0)

            stability_all(j,r) = -1;      % Saddle

        elseif all(real(lambda) > 0)

            stability_all(j,r) = 0;       % Unstable

        else

            stability_all(j,r) = NaN;     % Non-hyperbolic

        end

    end

end


%% Plot all equilibrium points with stability

figure;
hold on;

I_plot = repmat(I_values(:),1,max_eq);

stable = stability_all == 1;
saddle = stability_all == -1;
unstable = stability_all == 0;

h1 = plot(I_plot(stable), Veq_all(stable), ...
          'o','MarkerSize',5,'LineWidth',1.5);

h2 = plot(I_plot(saddle), Veq_all(saddle), ...
          'x','MarkerSize',7,'LineWidth',1.5);

h3 = plot(I_plot(unstable), Veq_all(unstable), ...
          's','MarkerSize',5,'LineWidth',1.5);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('V_{eq} (mV)');
title('Equilibrium Points and Stability');

grid on;

legend([h1 h2 h3], ...
       {'Stable','Saddle','Unstable'}, ...
       'Location','best');

hold off;


%% Plot equilibrium w versus current

figure;
hold on;

plot(I_values,Veq_all(:,1),'LineWidth',2);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('w_{eq}');

title('Equilibrium w versus Applied Current');

grid on;

hold off;


fprintf('Question 16\n');
%% Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;
EL  = -49.4;

phi = 1;

Iext = 10;

%% Rate functions

alpha_n = @(V) ...
    0.01*(V+50)./(1-exp(-(V+50)/10));

beta_n = @(V) ...
    0.125*exp(-(V+60)/80);

alpha_m = @(V) ...
    0.1*(V+35)./(1-exp(-(V+35)/10));

beta_m = @(V) ...
    4*exp(-(V+60)/18);

alpha_h = @(V) ...
    0.07*exp(-(V+60)/20);

beta_h = @(V) ...
    1./(1+exp(-(V+30)/10));


%% Steady-state values of m and h

m_inf = @(V) alpha_m(V)./(alpha_m(V)+beta_m(V));

h_inf = @(V) alpha_h(V)./(alpha_h(V)+beta_h(V));


%% Initial conditions

V0 = -60;

n0 = alpha_n(V0)/(alpha_n(V0)+beta_n(V0));

Y0 = [V0; n0];


%% Reduced HH model

MLE = @(t,Y) [

    (Iext ...
    - GK*Y(2)^4*(Y(1)-EK) ...
    - GNa*m_inf(Y(1))^3*h_inf(Y(1))*(Y(1)-ENa) ...
    - GL*(Y(1)-EL))/C;

    phi*(alpha_n(Y(1))*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2))

];


%% Solve

tspan = [0 100];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[t,Y] = ode15s(MLE,tspan,Y0,options);


%% Extract variables

V = Y(:,1);
n = Y(:,2);


%% Plot membrane potential

figure;

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');

title('Reduced HH Model: Membrane Potential');

grid on;


%% Plot n

figure;

plot(t,n,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('n');

title('Reduced HH Model: Potassium Activation');

grid on;


fprintf('Question 17\n');
%% Parameters

C = 1;                  % uF/cm^2

GK  = 36;               % mS/cm^2
GNa = 120;              % mS/cm^2
GL  = 0.3;              % mS/cm^2

EK  = -72;              % mV
ENa = 55;               % mV

%% Temperature

T = 6.3;
Q = 3;

phi = Q^((T-6.3)/10);

%% Determine EL for Vrest = -60 mV

Vrest = -60;
Irest = 0;

% Steady-state gating variables

ninf = alpha_n_safe(Vrest,phi) / ...
       (alpha_n_safe(Vrest,phi) + beta_n(Vrest));

minf = alpha_m_safe(Vrest,phi) / ...
       (alpha_m_safe(Vrest,phi) + beta_m(Vrest));

hinf = alpha_h(Vrest) / ...
       (alpha_h(Vrest) + beta_h(Vrest));

% Calculate EL

EL = Vrest - ...
     (-GK*ninf^4*(Vrest-EK) ...
      -GNa*minf^3*hinf*(Vrest-ENa) ...
      + Irest)/GL;

fprintf('Temperature = %.2f C\n',T);
fprintf('phi = %.4f\n',phi);

fprintf('\nResting-state gating variables:\n');
fprintf('n_inf = %.6f\n',ninf);
fprintf('m_inf = %.6f\n',minf);
fprintf('h_inf = %.6f\n',hinf);

fprintf('\nCalculated leak reversal potential:\n');
fprintf('EL = %.6f mV\n',EL);


%% Now apply external current

Iext = @(t) -3*(t <= 20);

%% Initial condition

% Start exactly at the resting state

Y0 = [Vrest; ninf; minf; hinf];


%% Hodgkin-Huxley equations

HH = @(t,Y) [

    % dV/dt

    (Iext(t) ...
    - GK*Y(2)^4*(Y(1)-EK) ...
    - GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    - GL*(Y(1)-EL))/C;


    % dn/dt

    phi*( ...
    alpha_n_safe(Y(1),phi)*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2));


    % dm/dt

    phi*( ...
    alpha_m_safe(Y(1),phi)*(1-Y(3)) ...
    - beta_m(Y(1))*Y(3));


    % dh/dt

    phi*( ...
    alpha_h(Y(1))*(1-Y(4)) ...
    - beta_h(Y(1))*Y(4))
];


%% Simulation

tspan = [0 1000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);


%% Use stiff solver ode15s

[t,Y] = ode15s(HH,tspan,Y0,options);


%% Extract variables

V = Y(:,1);
n = Y(:,2);
m = Y(:,3);
h = Y(:,4);


%% Plot membrane potential

figure;

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');

title('Anode-Break Excitation: -3 \muA/cm^2 for 20 ms');

grid on;


%% Plot gating variables

figure;

plot(t,n,'LineWidth',1.5);
hold on;

plot(t,m,'LineWidth',1.5);
plot(t,h,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Gating variable');

title('HH Gating Variables');

legend('n','m','h');

grid on;

hold off;


%% Safe alpha_n function

function a = alpha_n_safe(V,phi)

    if abs(V + 50) < 1e-8

        % Limit as V -> -50
        a = 0.1*phi;

    else

        a = -0.01*phi*(V+50) / ...
            (exp(-(V+50)/10)-1);

    end

end


%% Safe alpha_m function

function a = alpha_m_safe(V,phi)

    if abs(V + 35) < 1e-8

        % Limit as V -> -35
        a = phi;

    else

        a = -0.1*phi*(V+35) / ...
            (exp(-(V+35)/10)-1);

    end

end


%% beta_n

function b = beta_n(V)

    b = 0.125*exp(-(V+60)/80);

end


%% beta_m

function b = beta_m(V)

    b = 4*exp(-(V+60)/18);

end


%% alpha_h

function a = alpha_h(V)

    a = 0.07*exp(-(V+60)/20);

end


%% beta_h

function b = beta_h(V)

    b = 1./(exp(-(V+30)/10)+1);

end



fprintf('Question 18\n');
%% HH Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;
EL  = -49.4;

phi = 1;

%% Rate functions

alpha_n = @(V) alpha_n_safe(V,phi);
beta_n  = @(V) 0.125*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V,phi);
beta_m  = @(V) 4*exp(-(V+60)/18);

alpha_h = @(V) 0.07*exp(-(V+60)/20);
beta_h  = @(V) 1./(exp(-(V+30)/10)+1);

%% Steady-state gating variables

ninf = @(V) alpha_n(V)./(alpha_n(V)+beta_n(V));
minf = @(V) alpha_m(V)./(alpha_m(V)+beta_m(V));
hinf = @(V) alpha_h(V)./(alpha_h(V)+beta_h(V));
%% Voltage range and m range

V = linspace(-100,60,500);
m = linspace(0,1,500);

[VV,MM] = meshgrid(V,m);

%% =========================================================
% CASE 1: n and h fixed at their resting values
% ==========================================================

Vrest = -60;

n_fixed1 = ninf(Vrest);
h_fixed1 = hinf(Vrest);
Vdot1 = ( ...
    -GK*n_fixed1^4.*(VV-EK) ...
    -GNa*MM.^3*h_fixed1.*(VV-ENa) ...
    -GL*(VV-EL) ) / C;

% m equation = 0
Mdot1 = alpha_m(VV).*(1-MM) - beta_m(VV).*MM;

figure;
hold on;

contour(VV,MM,Vdot1,[0 0],'LineWidth',2);
contour(VV,MM,Mdot1,[0 0],'LineWidth',2,'LineColor','r');

xlabel('V (mV)');
ylabel('m');
title('V-m Phase Plane: n,h at Rest');

legend('V-nullcline','m-nullcline');

grid on;
xlim([-100 60]);
ylim([0 1]);
%% Find equilibrium points

F1 = @(Veq) ...
    -GK*n_fixed1^4*(Veq-EK) ...
    -GNa*minf(Veq)^3*h_fixed1*(Veq-ENa) ...
    -GL*(Veq-EL);

roots1 = [];

Vscan = linspace(-100,60,3000);
Fscan = arrayfun(F1,Vscan);

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1) < 0

        r = fzero(F1,[Vscan(i),Vscan(i+1)]);

        if isempty(roots1) || all(abs(roots1-r)>1e-5)
            roots1(end+1) = r;
        end

    end
end

mroots1 = minf(roots1);
plot(roots1,mroots1,'ko','MarkerFaceColor','k','MarkerSize',8);

fprintf('\nCASE 1: n,h at rest\n');
fprintf('n_fixed = %.6f\n',n_fixed1);
fprintf('h_fixed = %.6f\n',h_fixed1);
fprintf('Equilibrium points:\n');

for i = 1:length(roots1)
    fprintf('V = %.6f mV, m = %.6f\n', ...
        roots1(i),mroots1(i));
end

hold off;
%% =========================================================
% CASE 2: n and h at end of anodal stimulus
% ==========================================================

% ---------------------------------------------------------
% Obtain n and h at t = 20 ms from your Q17 simulation
% ---------------------------------------------------------

Iext = @(t) -3*(t <= 20);

Y0 = [Vrest; ninf(Vrest); minf(Vrest); hinf(Vrest)];

HH = @(t,Y) [

    (Iext(t) ...
    - GK*Y(2)^4*(Y(1)-EK) ...
    - GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    - GL*(Y(1)-EL))/C;

    alpha_n(Y(1))*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2);

    alpha_m(Y(1))*(1-Y(3)) ...
    - beta_m(Y(1))*Y(3);

    alpha_h(Y(1))*(1-Y(4)) ...
    - beta_h(Y(1))*Y(4)
];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[t,Y] = ode15s(HH,[0 20],Y0,options);

% Values at the end of the anodal stimulus

n_fixed2 = Y(end,2);
h_fixed2 = Y(end,4);

fprintf('\nCASE 2: n,h at end of anodal stimulus\n');
fprintf('n_fixed = %.6f\n',n_fixed2);
fprintf('h_fixed = %.6f\n',h_fixed2);

%% Phase plane for Case 2

Vdot2 = ( ...
    -GK*n_fixed2^4.*(VV-EK) ...
    -GNa*MM.^3*h_fixed2.*(VV-ENa) ...
    -GL*(VV-EL) ) / C;

Mdot2 = alpha_m(VV).*(1-MM) - beta_m(VV).*MM;

figure;
hold on;

contour(VV,MM,Vdot2,[0 0],'LineWidth',2);
contour(VV,MM,Mdot2,[0 0],'LineWidth',2);

xlabel('V (mV)');
ylabel('m');

title('V-m Phase Plane: n,h at End of Anodal Stimulus');

legend('V-nullcline','m-nullcline');

grid on;
xlim([-100 60]);
ylim([0 1]);

%% Find equilibria for Case 2

F2 = @(Veq) ...
    -GK*n_fixed2^4*(Veq-EK) ...
    -GNa*minf(Veq)^3*h_fixed2*(Veq-ENa) ...
    -GL*(Veq-EL);

roots2 = [];

Fscan = arrayfun(F2,Vscan);

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1) < 0

        r = fzero(F2,[Vscan(i),Vscan(i+1)]);

        if isempty(roots2) || all(abs(roots2-r)>1e-5)
            roots2(end+1) = r;
        end

    end

end

mroots2 = minf(roots2);

plot(roots2,mroots2,'ko','MarkerFaceColor','k','MarkerSize',8);

fprintf('\nCase 2 equilibrium points:\n');

for i = 1:length(roots2)
    fprintf('V = %.6f mV, m = %.6f\n', ...
        roots2(i),mroots2(i));
end

hold off;
%% =========================================================
% Safe alpha functions
% ==========================================================

function a = alpha_n_safe(V,phi)

    a = -0.01*phi.*(V+50) ./ ...
        (exp(-(V+50)/10)-1);

    idx = abs(V+35) < 1e-8;

    a(idx) = 0.1* phi;

end


function a = alpha_m_safe(V,phi)

    a = -0.1*phi.*(V+35) ./ ...
        (exp(-(V+35)/10)-1);

    idx = abs(V+35) < 1e-8;

    a(idx) = phi;

end