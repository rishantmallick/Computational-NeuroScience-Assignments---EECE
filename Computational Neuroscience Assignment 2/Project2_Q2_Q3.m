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