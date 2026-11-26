clear;
clc;
close all;

%% ========================================================================
% QUESTION 2 AND QUESTION 3
% Morris-Lecar Nullclines, Direction Field, Equilibrium and Stability
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 2 AND QUESTION 3\n');
fprintf('============================================================\n');

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

%% Anonymous Functions
minf = @(V) 0.5*(1 + tanh((V-V1)/V2));
winf = @(V) 0.5*(1 + tanh((V-V3)/V4));

V_null = @(V) ...
    (-gca*minf(V).*(V-vca) ...
     -gl*(V-vl) + Iext) ...
     ./ (gk*(V-vk));

tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Equilibrium using fzero

F = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + Iext;

Veq = fzero(F,-60);
weq = winf(Veq);

fprintf('\nMorris-Lecar Equilibrium Point:\n');
fprintf('Veq = %.6f mV\n',Veq);
fprintf('weq = %.6f\n',weq);

%% Equilibrium using fsolve

f = @(x) [ ...
    (-gca*minf(x(1))*(x(1)-vca) ...
     -gk*x(2)*(x(1)-vk) ...
     -gl*(x(1)-vl) + Iext)/C;

    phi*(winf(x(1))-x(2))/tauw(x(1))
];

options = optimoptions('fsolve','Display','off');

eq = fsolve(f,[-60;0.01],options);

Veq = eq(1);
weq = eq(2);

fprintf('\nSolved with simultaneous equilibrium:\n');
fprintf('Veq = %.6f mV\n',Veq);
fprintf('weq = %.6f\n',weq);

%% Grid and vector calculations

V_curves = linspace(-100,100,500);

W_w = winf(V_curves);
W_V = V_null(V_curves);

[Vgrid,Wgrid] = meshgrid(-90:8:50,-0.1:0.08:1.1);

m = minf(Vgrid);
w_inf = winf(Vgrid);
tau = tauw(Vgrid);

dV = (-gca*m.*(Vgrid-vca) ...
      -gk*Wgrid.*(Vgrid-vk) ...
      -gl*(Vgrid-vl) + Iext)/C;

dw = phi*(w_inf-Wgrid)./tau;

L = sqrt(dV.^2 + dw.^2);
L(L == 0) = 1;

dV_plot = dV./L;
dw_plot = dw./L;

%% Figure 1: Nullclines

figure('Name','Question 2 - Nullclines');

plot(V_curves,W_V,'r-','LineWidth',2.5, ...
    'DisplayName','V-nullcline');
hold on;

plot(V_curves,W_w,'b-','LineWidth',2.5, ...
    'DisplayName','w-nullcline');

plot(Veq,weq,'ko','MarkerSize',9, ...
    'MarkerFaceColor','y','LineWidth',2, ...
    'DisplayName','Equilibrium');

xlabel('V (mV)');
ylabel('w');
title('Question 2: Morris-Lecar Nullclines');

legend('Location','best');
grid on;

xlim([-90 50]);
ylim([-0.1 1.1]);

%% Figure 2: Direction Field

figure('Name','Question 3 - Direction Field');

quiver(Vgrid,Wgrid,dV_plot,dw_plot,0.4, ...
    'k','LineWidth',1.0);

xlabel('V (mV)');
ylabel('w');
title('Question 3: Morris-Lecar Vector Flow Field');

grid on;

xlim([-90 50]);
ylim([-0.1 1.1]);

%% Figure 3: Combined Phase Plane

figure('Name','Question 3 - Combined Phase Plane');

quiver(Vgrid,Wgrid,dV_plot,dw_plot,0.4, ...
    'k','LineWidth',1.0, ...
    'HandleVisibility','off');

hold on;

plot(V_curves,W_V,'r-','LineWidth',2.5, ...
    'DisplayName','V-nullcline');

plot(V_curves,W_w,'b-','LineWidth',2.5, ...
    'DisplayName','w-nullcline');

plot(Veq,weq,'ko','MarkerSize',10, ...
    'MarkerFaceColor','y','LineWidth',2, ...
    'DisplayName',sprintf('Equilibrium: (%.2f, %.2f)',Veq,weq));

xlabel('V (mV)');
ylabel('w');
title('Question 3: Morris-Lecar Phase Plane');

legend('Location','best');
grid on;

xlim([-90 50]);
ylim([-0.1 1.1]);

%% Jacobian and Eigenvalues

syms V w

m_sym = 0.5*(1+tanh((V-V1)/V2));
winf_sym = 0.5*(1+tanh((V-V3)/V4));
tauw_sym = 1/cosh((V-V3)/(2*V4));

f1 = (-gca*m_sym*(V-vca) ...
      -gk*w*(V-vk) ...
      -gl*(V-vl) + Iext)/C;

f2 = phi*(winf_sym-w)/tauw_sym;

J = jacobian([f1;f2],[V w]);

J_eq = double(subs(J,[V w],[Veq weq]));

eigenvalues = eig(J_eq);

fprintf('\nJacobian at equilibrium:\n');
disp(J_eq);

fprintf('Eigenvalues:\n');
disp(eigenvalues);


%% ========================================================================
% QUESTION 6
% Initial Depolarization and Threshold-like Response
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 6\n');
fprintf('============================================================\n');

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
Iext = 0;

%% Functions

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Resting equilibrium

F = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + Iext;

Veq = fzero(F,-60);
weq = winf(Veq);

fprintf('\nResting equilibrium:\n');
fprintf('Veq = %.6f mV\n',Veq);
fprintf('weq = %.6f\n',weq);

%% Morris-Lecar ODE

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;

    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

%% Initial depolarizations

V0_values = Veq + ...
    [1 2 3 4 5 6 7 8 9 10 12 15 20];

tspan = [0 300];

Vmax_values = zeros(size(V0_values));

trajectories = cell(length(V0_values),1);

%% Simulate

for i = 1:length(V0_values)

    V0 = V0_values(i);

    Y0 = [V0;weq];

    [t,Y] = ode45(MLE,tspan,Y0);

    V = Y(:,1);
    w = Y(:,2);

    Vmax_values(i) = max(V);

    trajectories{i}.V = V;
    trajectories{i}.w = w;
    trajectories{i}.t = t;

    fprintf('V0 = %.6f mV   Vmax = %.6f mV\n', ...
        V0,Vmax_values(i));

end

%% Phase-plane trajectories

figure('Name','Question 6 - Phase Plane');

hold on;

for i = 1:length(V0_values)

    plot(trajectories{i}.V, ...
         trajectories{i}.w, ...
         'LineWidth',1.2, ...
         'DisplayName',sprintf('V_0 = %.3f',V0_values(i)));

end

plot(Veq,weq,'ko','MarkerSize',8, ...
    'MarkerFaceColor','y','LineWidth',2, ...
    'DisplayName','Equilibrium');

xlabel('V (mV)');
ylabel('w');

title('Question 6: Morris-Lecar Phase-Plane Trajectories');

legend('Location','best');
grid on;

%% Maximum amplitude

figure('Name','Question 6 - Maximum Amplitude');

plot(V0_values,Vmax_values,'ko-', ...
    'LineWidth',1.5,'MarkerSize',5);

xlabel('Initial voltage V_0 (mV)');
ylabel('Maximum voltage V_{max} (mV)');

title('Question 6: Maximum Action-Potential Amplitude');

grid on;

%% Finer threshold search

V0_values = Veq + (0:60);

Vthreshold_values = zeros(size(V0_values));

trajectories = cell(length(V0_values),1);

for i = 1:length(V0_values)

    V0 = V0_values(i);

    Y0 = [V0;weq];

    [t,Y] = ode45(MLE,tspan,Y0);

    V = Y(:,1);
    w = Y(:,2);

    Vthreshold_values(i) = max(V);

    trajectories{i}.V = V;
    trajectories{i}.w = w;
    trajectories{i}.t = t;

    fprintf('V0 = %.6f mV   Vmax = %.6f mV\n', ...
        V0,Vthreshold_values(i));

end

figure('Name','Question 6 - Threshold Search');

plot(V0_values,Vthreshold_values,'ko-', ...
    'LineWidth',1.5,'MarkerSize',5);

xlabel('Initial voltage V_0 (mV)');
ylabel('Maximum voltage V_{max} (mV)');

title('Question 6: Maximum Potential vs Initial Voltage');

grid on;

%% Below and above threshold

V0_below = Veq + 40;
V0_above = Veq + 50;

figure('Name','Question 6 - Below Above Threshold');

hold on;

[t1,Y1] = ode45(MLE,tspan,[V0_below;weq]);
[t2,Y2] = ode45(MLE,tspan,[V0_above;weq]);

plot(Y1(:,1),Y1(:,2),'LineWidth',2, ...
    'DisplayName',sprintf('Below: V_0 = %.3f',V0_below));

plot(Y2(:,1),Y2(:,2),'LineWidth',2, ...
    'DisplayName',sprintf('Above: V_0 = %.3f',V0_above));

plot(Veq,weq,'ko','MarkerSize',8, ...
    'MarkerFaceColor','y','LineWidth',2, ...
    'DisplayName','Equilibrium');

xlabel('V (mV)');
ylabel('w');

title('Question 6: Phase Plane Below and Above Threshold');

legend('Location','best');
grid on;


%% ========================================================================
% QUESTION 5
% Effect of phi
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 5\n');
fprintf('============================================================\n');

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

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

F = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + Iext;

Veq = fzero(F,-60);
weq = winf(Veq);

tspan = [0 300];

Y0 = [Veq+1;weq];

%% phi = 0.01

phi = 0.01;

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;

    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

[t01,Y01] = ode45(MLE,tspan,Y0);

%% phi = 0.02

phi = 0.02;

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;

    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

[t02,Y02] = ode45(MLE,tspan,Y0);

%% phi = 0.04

phi = 0.04;

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;

    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

[t04,Y04] = ode45(MLE,tspan,Y0);

%% Action potential comparison

figure('Name','Question 5 - Action Potentials');

plot(t01,Y01(:,1),'LineWidth',1.5);
hold on;

plot(t02,Y02(:,1),'LineWidth',1.5);
plot(t04,Y04(:,1),'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');

title('Question 5: Action Potential for Different \phi');

legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

%% Phase plane comparison

figure('Name','Question 5 - Phase Plane');

plot(Y01(:,1),Y01(:,2),'LineWidth',1.5);
hold on;

plot(Y02(:,1),Y02(:,2),'LineWidth',1.5);
plot(Y04(:,1),Y04(:,2),'LineWidth',1.5);

plot(Veq,weq,'ko','MarkerSize',8, ...
    'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');

title('Question 5: Phase Plane for Different \phi');

legend('\phi = 0.01','\phi = 0.02', ...
       '\phi = 0.04','Equilibrium');

grid on;


%% ========================================================================
% QUESTION 7
% Different Initial Conditions at Iext = 86
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 7\n');
fprintf('============================================================\n');

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

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Equilibrium for Iext = 0

I0 = 0;

F0 = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + I0;

Veq0 = fzero(F0,-60);
weq0 = winf(Veq0);

fprintf('\nEquilibrium for Iext = 0:\n');
fprintf('Veq0 = %.6f mV\n',Veq0);
fprintf('weq0 = %.6f\n',weq0);

%% Equilibrium for Iext = 86

F86 = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + Iext;

Veq86 = fzero(F86,-30);
weq86 = winf(Veq86);

fprintf('\nEquilibrium for Iext = 86:\n');
fprintf('Veq86 = %.6f mV\n',Veq86);
fprintf('weq86 = %.6f\n',weq86);

%% Morris-Lecar equations

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;

    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

tspan = [0 1000];

%% Initial conditions

Y01 = [Veq0;weq0];
Y02 = [Veq86;weq86];
Y03 = [-27.9;0.17];

%% Simulate

[t1,Y1] = ode45(MLE,tspan,Y01);
[t2,Y2] = ode45(MLE,tspan,Y02);
[t3,Y3] = ode45(MLE,tspan,Y03);

%% Phase plane

figure('Name','Question 7 - Phase Plane');

hold on;

plot(Y1(:,1),Y1(:,2),'LineWidth',1.5);
plot(Y2(:,1),Y2(:,2),'LineWidth',1.5);
plot(Y3(:,1),Y3(:,2),'LineWidth',1.5);

plot(Veq86,weq86,'ko','MarkerSize',8, ...
    'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');

title('Question 7: Morris-Lecar Phase Plane, I_{ext}=86');

legend( ...
    'Trial 1: I_{ext}=0 equilibrium', ...
    'Trial 2: I_{ext}=86 equilibrium', ...
    'Trial 3: (-27.9, 0.17)', ...
    'Equilibrium for I_{ext}=86', ...
    'Location','best');

grid on;


%% ========================================================================
% QUESTION 8
% Backward Integration and Unstable Periodic Orbit
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 8\n');
fprintf('============================================================\n');

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

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

MLE = @(t,Y) [ ...
    (-gca*minf(Y(1))*(Y(1)-vca) ...
     -gk*Y(2)*(Y(1)-vk) ...
     -gl*(Y(1)-vl) + Iext)/C;

    phi*(winf(Y(1))-Y(2))/tauw(Y(1))
];

%% Equilibrium

F = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + Iext;

Veq = fzero(F,-30);
weq = winf(Veq);

fprintf('\nEquilibrium:\n');
fprintf('Veq = %.6f mV\n',Veq);
fprintf('weq = %.6f\n',weq);

%% Backward integration

Y0 = [-27.9;0.17];

tspan_backward = [0 -5000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

[tb,Yb] = ode45(MLE,tspan_backward,Y0,options);

%% Backward trajectory

figure('Name','Question 8 - Backward Trajectory');

plot(Yb(:,1),Yb(:,2),'LineWidth',1.2);
hold on;

plot(Veq,weq,'ko','MarkerSize',8, ...
    'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');

title('Question 8: Backward-Time Trajectory');

legend('Backward trajectory','Equilibrium');
grid on;

%% UPO

figure('Name','Question 8 - UPO');

plot(Yb(:,1),Yb(:,2),'LineWidth',1.5);
hold on;

plot(Veq,weq,'ko','MarkerSize',8, ...
    'MarkerFaceColor','y');

xlabel('V (mV)');
ylabel('w');

title('Question 8: Unstable Periodic Orbit');

legend('UPO','Equilibrium');
grid on;


%% ========================================================================
% QUESTION 9
% Firing Rate versus Applied Current
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 9\n');
fprintf('============================================================\n');

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

    Iext = I(i);

    %% MLE

    MLE = @(t,Y) [ ...
        (-gca*minf(Y(1))*(Y(1)-vca) ...
         -gk*Y(2)*(Y(1)-vk) ...
         -gl*(Y(1)-vl) + Iext)/C;

        phi*(winf(Y(1))-Y(2))/tauw(Y(1))
    ];

    %% Equilibrium

    F = @(V) ...
        -gca*minf(V).*(V-vca) ...
        -gk*winf(V).*(V-vk) ...
        -gl*(V-vl) + Iext;

    Veq = fzero(F,-30);
    weq = winf(Veq);

    %% Jacobian

    syms V w

    f1 = (-gca*minf(V)*(V-vca) ...
          -gk*w*(V-vk) ...
          -gl*(V-vl) + Iext)/C;

    f2 = phi*(winf(V)-w)/tauw(V);

    J = jacobian([f1;f2],[V w]);

    J_eq = double(subs(J,[V w],[Veq weq]));

    eigenvalues = eig(J_eq);

    fprintf('\nIext = %g\n',Iext);
    fprintf('Jacobian:\n');
    disp(J_eq);

    fprintf('Eigenvalues:\n');
    disp(eigenvalues);

    %% Simulation

    tspan = [0 3000];

    Y0 = [Veq;weq];

    options = odeset( ...
        'RelTol',1e-8, ...
        'AbsTol',1e-10);

    [t,Y] = ode45(MLE,tspan,Y0,options);

    Vtrace = Y(:,1);

    %% Voltage plot

    figure;

    plot(t,Vtrace);

    xlabel('Time (ms)');
    ylabel('Voltage (mV)');

    title(['Question 9: V vs t, Iext = ',num2str(Iext)]);

    grid on;

    filename = ['V_vs_t_Iext_',num2str(Iext),'_Q9.png'];

    exportgraphics(gcf,filename,'Resolution',300);

    %% Spike detection

    [pks,locs] = findpeaks(Vtrace,t, ...
        'MinPeakHeight',-50);

    %% ISI

    ISI = diff(locs);

    %% Firing rate

    if ~isempty(ISI)

        firing_rate = 1000/mean(ISI);

    else

        firing_rate = 0;

    end

    firingrate(i) = firing_rate;

    fprintf('Spike times (ms):\n');
    disp(locs);

    fprintf('ISI (ms):\n');
    disp(ISI);

    fprintf('Firing rate = %.4f Hz\n',firing_rate);

end

%% Firing-rate curve

figure('Name','Question 9 - Firing Rate');

plot(I,firingrate,'o-');

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Firing Rate (Hz)');

title('Question 9: Firing Rate vs Applied Current');

grid on;

exportgraphics(gcf, ...
    'Firing_Rate_vs_Applied_Current_Q9.png', ...
    'Resolution',300);


%% ========================================================================
% QUESTION 10
% Stable and Unstable Manifolds
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 10\n');
fprintf('============================================================\n');

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

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Equilibrium points

F = @(V) ...
    -gca*minf(V).*(V-vca) ...
    -gk*winf(V).*(V-vk) ...
    -gl*(V-vl) + Iext;

Vscan = linspace(-100,100,2000);

Fscan = F(Vscan);

rootsV = [];

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1) < 0

        root = fzero(F,[Vscan(i),Vscan(i+1)]);

        if isempty(rootsV) || ...
                all(abs(rootsV-root)>1e-5)

            rootsV(end+1) = root;

        end

    end

end

Veq = rootsV;
weq = winf(Veq);

disp('Equilibrium points:');
disp([Veq(:),weq(:)]);

%% Symbolic Jacobian

syms Vsym wsym

msym = 0.5*(1+tanh((Vsym-V1)/V2));
winfsym = 0.5*(1+tanh((Vsym-V3)/V4));
tauwsym = 1/cosh((Vsym-V3)/(2*V4));

f1 = (-gca*msym*(Vsym-vca) ...
      -gk*wsym*(Vsym-vk) ...
      -gl*(Vsym-vl) + Iext)/C;

f2 = phi*(winfsym-wsym)/tauwsym;

J = jacobian([f1;f2],[Vsym wsym]);

%% Stability

for k = 1:length(Veq)

    J_eq = double(subs(J, ...
        [Vsym wsym], ...
        [Veq(k) weq(k)]));

    lambda = eig(J_eq);

    fprintf('\nEquilibrium %d:\n',k);
    fprintf('Veq = %.6f, weq = %.6f\n', ...
        Veq(k),weq(k));

    disp('Jacobian:');
    disp(J_eq);

    disp('Eigenvalues:');
    disp(lambda);

    if all(real(lambda)<0)

        disp('Stable');

    elseif any(real(lambda)>0) && ...
           any(real(lambda)<0)

        disp('Saddle point');

    elseif all(real(lambda)>0)

        disp('Unstable');

    else

        disp('Non-hyperbolic');

    end

end

%% Nullclines

V = linspace(-100,100,2000);

V_null = ...
    (-gca*minf(V).*(V-vca) ...
     -gl*(V-vl) + Iext) ...
     ./ (gk*(V-vk));

w_null = winf(V);

%% Find saddle

saddle_index = [];

for k = 1:length(Veq)

    J_eq = double(subs(J, ...
        [Vsym wsym], ...
        [Veq(k) weq(k)]));

    lambda = eig(J_eq);

    if real(lambda(1))*real(lambda(2)) < 0

        saddle_index = k;

        break;

    end

end

Vs = Veq(saddle_index);
Ws = weq(saddle_index);

J_saddle = double(subs(J, ...
    [Vsym wsym], ...
    [Vs Ws]));

[VEC,D] = eig(J_saddle);

lambda = diag(D);

unstable_index = find(real(lambda)>0,1);
stable_index   = find(real(lambda)<0,1);

v_unstable = real(VEC(:,unstable_index));
v_stable   = real(VEC(:,stable_index));

%% ODE

odefun = @(t,y) [ ...
    (-gca*minf(y(1))*(y(1)-vca) ...
     -gk*y(2)*(y(1)-vk) ...
     -gl*(y(1)-vl) + Iext)/C;

    phi*(winf(y(1))-y(2))/tauw(y(1))
];

epsilon = 1e-3;

y0_u1 = [Vs;Ws] + epsilon*v_unstable;
y0_u2 = [Vs;Ws] - epsilon*v_unstable;

[t_u1,Y_u1] = ode45(odefun,[0 50],y0_u1);
[t_u2,Y_u2] = ode45(odefun,[0 50],y0_u2);

y0_s1 = [Vs;Ws] + epsilon*v_stable;
y0_s2 = [Vs;Ws] - epsilon*v_stable;

[t_s1,Y_s1] = ode45(odefun,[0 -50],y0_s1);
[t_s2,Y_s2] = ode45(odefun,[0 -50],y0_s2);

%% Nullclines and equilibria

figure('Name','Question 10 - Nullclines');

hold on;

plot(V,V_null,'LineWidth',2);
plot(V,w_null,'LineWidth',2);

plot(Veq,weq,'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

xlabel('V (mV)');
ylabel('w');

title('Question 10: Morris-Lecar Nullclines and Equilibria');

legend('V-nullcline', ...
       'w-nullcline', ...
       'Equilibrium points');

grid on;

xlim([-50 15]);
ylim([0 0.4]);

hold off;

%% Stable and unstable manifolds

figure('Name','Question 10 - Manifolds');

hold on;

hU = plot(Y_u1(:,1),Y_u1(:,2),'-', ...
    'LineWidth',2.5);

plot(Y_u2(:,1),Y_u2(:,2),'-', ...
    'LineWidth',2.5);

hS = plot(Y_s1(:,1),Y_s1(:,2),'--', ...
    'LineWidth',2.5);

plot(Y_s2(:,1),Y_s2(:,2),'--', ...
    'LineWidth',2.5);

hE = plot(Vs,Ws,'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

xlabel('V (mV)');
ylabel('w');

title('Question 10: Stable and Unstable Manifolds');

legend([hU hS hE], ...
       'Unstable manifold', ...
       'Stable manifold', ...
       'Saddle point');

grid on;

xlim([-50 15]);
ylim([0 0.4]);

hold off;


%% ========================================================================
% QUESTION 11
% Equilibrium Bifurcation and Firing Rate
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 11\n');
fprintf('============================================================\n');

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

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

%% Equilibrium versus current

Vscan = linspace(-100,100,2000);

I_values = linspace(30,50,200);

Veq_all = NaN(length(I_values),3);
weq_all = NaN(length(I_values),3);

for j = 1:length(I_values)

    Iext = I_values(j);

    F = @(V) ...
        -gca*minf(V).*(V-vca) ...
        -gk*winf(V).*(V-vk) ...
        -gl*(V-vl) + Iext;

    Fscan = F(Vscan);

    rootsV = [];

    for i = 1:length(Vscan)-1

        if Fscan(i)*Fscan(i+1)<0

            root = fzero(F,[Vscan(i),Vscan(i+1)]);

            if isempty(rootsV) || ...
                    all(abs(rootsV-root)>1e-5)

                rootsV(end+1) = root;

            end

        end

    end

    nroots = min(length(rootsV),3);

    if nroots > 0

        Veq_all(j,1:nroots) = rootsV(1:nroots);
        weq_all(j,1:nroots) = winf(rootsV(1:nroots));

    end

end

%% Narrow range 39 to 40

Vscan1 = linspace(-100,100,2000);

I_values1 = linspace(39,40,100);

Veq_all1 = NaN(length(I_values1),3);

weq_all1 = NaN(length(I_values1),3);

for j = 1:length(I_values1)

    Iext = I_values1(j);

    F = @(V) ...
        -gca*minf(V).*(V-vca) ...
        -gk*winf(V).*(V-vk) ...
        -gl*(V-vl) + Iext;

    Fscan = F(Vscan1);

    rootsV1 = [];

    for i = 1:length(Vscan1)-1

        if Fscan(i)*Fscan(i+1)<0

            root = fzero(F,[Vscan1(i),Vscan1(i+1)]);

            if isempty(rootsV1) || ...
                    all(abs(rootsV1-root)>1e-5)

                rootsV1(end+1) = root;

            end

        end

    end

    nroots = min(length(rootsV1),3);

    if nroots > 0

        Veq_all1(j,1:nroots) = rootsV1(1:nroots);
        weq_all1(j,1:nroots) = winf(rootsV1(1:nroots));

    end

end

%% Plot equilibria

figure('Name','Question 11 - Equilibrium Voltage');

hold on;

plot(I_values1,Veq_all1(:,1),'LineWidth',2);
plot(I_values1,Veq_all1(:,2),'LineWidth',2);
plot(I_values1,Veq_all1(:,3),'LineWidth',2);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('V_{eq} (mV)');

title('Question 11: Equilibrium Points for I_{ext}=39 to 40');

legend('Equilibrium 1', ...
       'Equilibrium 2', ...
       'Equilibrium 3');

grid on;

xlim([39 40]);

hold off;

%% Firing rate

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

I_values = 30:0.25:45;

firing_rate = zeros(size(I_values));

tspan = [0 2000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

for j = 1:length(I_values)

    Iext = I_values(j);

    MLE = @(t,y) [ ...
        (-gca*minf(y(1))*(y(1)-vca) ...
         -gk*y(2)*(y(1)-vk) ...
         -gl*(y(1)-vl) + Iext)/C;

        phi*(winf(y(1))-y(2))/tauw(y(1))
    ];

    Y0 = [-41.8452;0.002047];

    [t,Y] = ode45(MLE,tspan,Y0,options);

    Vtrace = Y(:,1);

    [pks,locs] = findpeaks(Vtrace,t, ...
        'MinPeakHeight',-50);

    ISI = diff(locs);

    if ~isempty(ISI)

        firing_rate(j) = 1000/mean(ISI);

    else

        firing_rate(j) = 0;

    end

end

disp('Iext       Firing Rate (Hz)');
disp([I_values(:),firing_rate(:)]);

figure('Name','Question 11 - Firing Rate');

plot(I_values,firing_rate,'o-', ...
    'LineWidth',2);

xlabel('Applied current I_{ext} (\muA/cm^2)');
ylabel('Firing rate (Hz)');

title('Question 11: Firing Rate vs Applied Current');

grid on;


%% ========================================================================
% QUESTION 12
% Classical Hodgkin-Huxley Model
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 12\n');
fprintf('============================================================\n');

%% Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;

T = 6.3;
Q = 3;

phi = Q^((T-6.3)/10);

Iext = 10;

EL = -49.4;

%% Rate functions

alpha_n = @(V) alpha_n_safe(V);
beta_n  = @(V) 0.125*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V);
beta_m  = @(V) 4*exp(-(V+60)/18);

alpha_h = @(V) 0.07*exp(-(V+60)/20);
beta_h  = @(V) 1./(exp(-(V+30)/10)+1);

%% Initial voltage

V0 = -60;

an0 = alpha_n(V0);
am0 = alpha_m(V0);

ah0 = alpha_h(V0);

bn0 = beta_n(V0);
bm0 = beta_m(V0);
bh0 = beta_h(V0);

n0 = an0/(an0+bn0);
m0 = am0/(am0+bm0);
h0 = ah0/(ah0+bh0);

Y0 = [V0;n0;m0;h0];

%% HH equations

HH = @(t,Y) [ ...
    (Iext ...
    -GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    -GL*(Y(1)-EL))/C;

    phi*( ...
    alpha_n(Y(1))*(1-Y(2)) ...
    -beta_n(Y(1))*Y(2));

    phi*( ...
    alpha_m(Y(1))*(1-Y(3)) ...
    -beta_m(Y(1))*Y(3));

    phi*( ...
    alpha_h(Y(1))*(1-Y(4)) ...
    -beta_h(Y(1))*Y(4))
];

tspan = [0 1000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

[t,Y] = ode45(HH,tspan,Y0,options);

V = Y(:,1);
n = Y(:,2);
m = Y(:,3);
h = Y(:,4);

GK_t = GK*n.^4;
GNa_t = GNa*m.^3.*h;

%% Voltage

figure('Name','Question 12 - Voltage');

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');

title(['Question 12: HH Model, I_{ext}=', ...
       num2str(Iext),' \muA/cm^2']);

grid on;

%% Gating variables

figure('Name','Question 12 - Gating');

plot(t,n,'LineWidth',1.5);
hold on;

plot(t,m,'LineWidth',1.5);
plot(t,h,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Gating variables');

title('Question 12: HH Gating Variables');

legend('n','m','h');

grid on;

hold off;

%% Conductances

figure('Name','Question 12 - Conductances');

plot(t,GK_t,'LineWidth',1.5);
hold on;

plot(t,GNa_t,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Conductance (mS/cm^2)');

title('Question 12: HH Conductances');

legend('G_K','G_{Na}');

grid on;

hold off;

fprintf('\nTemperature = %.1f degrees C\n',T);
fprintf('Phi = %.4f\n',phi);

fprintf('\nInitial conditions:\n');
fprintf('V0 = %.4f mV\n',V0);
fprintf('n0 = %.4f\n',n0);
fprintf('m0 = %.4f\n',m0);
fprintf('h0 = %.4f\n',h0);


%% ========================================================================
% QUESTION 13
% HH Model with EL chosen for Vrest = -60 mV
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 13\n');
fprintf('============================================================\n');

%% Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;

T = 6.3;
Q = 3;

phi = Q^((T-6.3)/10);

%% Resting potential

Vrest = -60;
Irest = 0;

ninf = alpha_n_safe(Vrest) / ...
       (alpha_n_safe(Vrest)+beta_n(Vrest));

minf = alpha_m_safe(Vrest) / ...
       (alpha_m_safe(Vrest)+beta_m(Vrest));

hinf = alpha_h(Vrest) / ...
       (alpha_h(Vrest)+beta_h(Vrest));

%% Calculate EL

EL = Vrest - ...
     (-GK*ninf^4*(Vrest-EK) ...
      -GNa*minf^3*hinf*(Vrest-ENa) ...
      + Irest)/GL;

fprintf('\nTemperature = %.2f C\n',T);
fprintf('phi = %.4f\n',phi);

fprintf('\nResting-state gating variables:\n');
fprintf('n_inf = %.6f\n',ninf);
fprintf('m_inf = %.6f\n',minf);
fprintf('h_inf = %.6f\n',hinf);

fprintf('\nCalculated leak reversal potential:\n');
fprintf('EL = %.6f mV\n',EL);

%% External current

Iext = 10;

Y0 = [Vrest;ninf;minf;hinf];

%% HH equations

HH = @(t,Y) [ ...
    (Iext ...
    -GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    -GL*(Y(1)-EL))/C;

    phi*( ...
    alpha_n_safe(Y(1))*(1-Y(2)) ...
    -beta_n(Y(1))*Y(2));

    phi*( ...
    alpha_m_safe(Y(1))*(1-Y(3)) ...
    -beta_m(Y(1))*Y(3));

    phi*( ...
    alpha_h(Y(1))*(1-Y(4)) ...
    -beta_h(Y(1))*Y(4))
];

tspan = [0 1000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

[t,Y] = ode15s(HH,tspan,Y0,options);

V = Y(:,1);
n = Y(:,2);
m = Y(:,3);
h = Y(:,4);

%% Voltage

figure('Name','Question 13 - Voltage');

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');

title(['Question 13: HH Model, I_{ext}=', ...
       num2str(Iext),' \muA/cm^2']);

grid on;

%% Gating

figure('Name','Question 13 - Gating');

plot(t,n,'LineWidth',1.5);
hold on;

plot(t,m,'LineWidth',1.5);
plot(t,h,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Gating variable');

title('Question 13: HH Gating Variables');

legend('n','m','h');

grid on;

hold off;


%% ========================================================================
% QUESTION 14
% HH Resting Equilibrium, Stability and Threshold
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 14\n');
fprintf('============================================================\n');

%% Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;

T = 6.3;
Q = 3;

EL = -49.4;

phi = Q^((T-6.3)/10);

%% Rate functions

alpha_n = @(V) alpha_n_safe(V);
beta_n = @(V) beta_n(V);

alpha_m = @(V) alpha_m_safe(V);
beta_m = @(V) beta_m(V);

alpha_h = @(V) 0.07*exp(-(V+60)/20);
beta_h = @(V) 1./(exp(-(V+30)/10)+1);

%% Steady-state gating variables

n_inf = @(V) alpha_n(V)./(alpha_n(V)+beta_n(V));
m_inf = @(V) alpha_m(V)./(alpha_m(V)+beta_m(V));
h_inf = @(V) alpha_h(V)./(alpha_h(V)+beta_h(V));

%% Resting equilibrium

Iext = 0;

F = @(V) ...
    -GK*n_inf(V).^4.*(V-EK) ...
    -GNa*m_inf(V).^3.*h_inf(V).*(V-ENa) ...
    -GL*(V-EL) + Iext;

Vscan = linspace(-100,100,3000);

Fscan = F(Vscan);

rootsV = [];

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1)<0

        root = fzero(F,[Vscan(i),Vscan(i+1)]);

        if isempty(rootsV) || ...
                all(abs(rootsV-root)>1e-5)

            rootsV(end+1) = root;

        end

    end

end

Vrest = rootsV(1);

nrest = n_inf(Vrest);
mrest = m_inf(Vrest);
hrest = h_inf(Vrest);

fprintf('\nRESTING EQUILIBRIUM\n');
fprintf('Vrest = %.6f mV\n',Vrest);
fprintf('nrest = %.6f\n',nrest);
fprintf('mrest = %.6f\n',mrest);
fprintf('hrest = %.6f\n',hrest);

%% Jacobian

syms Vs ns ms hs

f1 = (-GK*ns^4*(Vs-EK) ...
      -GNa*ms^3*hs*(Vs-ENa) ...
      -GL*(Vs-EL) + Iext)/C;

alpha_n_sym = ...
    -0.01*(Vs+50)/(exp(-(Vs+50)/10)-1);

beta_n_sym = ...
    0.125*exp(-(Vs+60)/80);

alpha_m_sym = ...
    -0.1*(Vs+35)/(exp(-(Vs+35)/10)-1);

beta_m_sym = ...
    4*exp(-(Vs+60)/18);

alpha_h_sym = ...
    0.07*exp(-(Vs+60)/20);

beta_h_sym = ...
    1/(exp(-(Vs+30)/10)+1);

f2 = alpha_n_sym*(1-ns) ...
     -beta_n_sym*ns;

f3 = alpha_m_sym*(1-ms) ...
     -beta_m_sym*ms;

f4 = alpha_h_sym*(1-hs) ...
     -beta_h_sym*hs;

J = jacobian([f1;f2;f3;f4], ...
             [Vs ns ms hs]);

Jrest = double(subs(J, ...
    [Vs ns ms hs], ...
    [Vrest nrest mrest hrest]));

eigenvalues = eig(Jrest);

fprintf('\nSTABILITY AT REST\n');

disp('Jacobian:');
disp(Jrest);

disp('Eigenvalues:');
disp(eigenvalues);

if all(real(eigenvalues)<0)

    fprintf('Resting equilibrium is STABLE.\n');

elseif any(real(eigenvalues)>0)

    fprintf('Resting equilibrium is UNSTABLE.\n');

else

    fprintf('Resting equilibrium is NON-HYPERBOLIC.\n');

end

%% HH equations

HH = @(t,Y) [ ...
    (-GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    -GL*(Y(1)-EL) + Iext)/C;

    alpha_n(Y(1))*(1-Y(2)) ...
    -beta_n(Y(1))*Y(2);

    alpha_m(Y(1))*(1-Y(3)) ...
    -beta_m(Y(1))*Y(3);

    alpha_h(Y(1))*(1-Y(4)) ...
    -beta_h(Y(1))*Y(4)
];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

%% Threshold search

deltaV_values = 0:0.1:40;

threshold_found = false;

spike_result = zeros(size(deltaV_values));

for k = 1:length(deltaV_values)

    deltaV = deltaV_values(k);

    Y0 = [ ...
        Vrest + deltaV;
        nrest;
        mrest;
        hrest
    ];

    [t,Y] = ode15s(HH,[0 100],Y0,options);

    V = Y(:,1);

    if max(V)>0

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

    fprintf('\nThreshold depolarization = %.2f mV\n', ...
        threshold_deltaV);

    fprintf('Threshold voltage = %.2f mV\n', ...
        threshold_voltage);

else

    fprintf('\nNo spike found in tested range.\n');

end

%% Threshold plot

figure('Name','Question 14 - Threshold');

plot(deltaV_values,spike_result,'o-', ...
    'LineWidth',1.5);

xlabel('\DeltaV (mV)');
ylabel('Spike');

title('Question 14: Threshold Depolarization');

grid on;

%% Voltage at threshold

if threshold_found

    figure('Name','Question 14 - Threshold Voltage');

    plot(t_threshold,V_threshold, ...
        'LineWidth',1.5);

    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');

    title(['Question 14: Response at Threshold, \DeltaV = ', ...
        num2str(threshold_deltaV),' mV']);

    grid on;

end


%% ========================================================================
% QUESTION 15
% Morris-Lecar Equilibrium Stability for Iext = 8 to 12
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 15\n');
fprintf('============================================================\n');

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

minf = @(V) 0.5*(1+tanh((V-V1)/V2));
winf = @(V) 0.5*(1+tanh((V-V3)/V4));
tauw = @(V) 1./cosh((V-V3)/(2*V4));

I_values = linspace(8,12,200);

Vscan = linspace(-100,100,3000);

max_eq = 3;

Veq_all = NaN(length(I_values),max_eq);
weq_all = NaN(length(I_values),max_eq);

stability_all = NaN(length(I_values),max_eq);

%% Symbolic Jacobian

syms V w

m = 0.5*(1+tanh((V-V1)/V2));

w_inf = 0.5*(1+tanh((V-V3)/V4));

tau = 1/cosh((V-V3)/(2*V4));

f1 = (-gca*m*(V-vca) ...
      -gk*w*(V-vk) ...
      -gl*(V-vl))/C;

f2 = phi*(w_inf-w)/tau;

J = jacobian([f1;f2],[V w]);

%% Find equilibria

for j = 1:length(I_values)

    Iext = I_values(j);

    F = @(V) ...
        -gca*minf(V).*(V-vca) ...
        -gk*winf(V).*(V-vk) ...
        -gl*(V-vl) + Iext;

    Fscan = F(Vscan);

    rootsV = [];

    for i = 1:length(Vscan)-1

        if Fscan(i)*Fscan(i+1)<0

            root = fzero(F,[Vscan(i),Vscan(i+1)]);

            if isempty(rootsV) || ...
                    all(abs(rootsV-root)>1e-5)

                rootsV(end+1) = root;

            end

        end

    end

    nroots = min(length(rootsV),max_eq);

    for r = 1:nroots

        Veq = rootsV(r);
        weq = winf(Veq);

        J_eq = double(subs(J,[V w],[Veq weq]));

        lambda = eig(J_eq);

        Veq_all(j,r) = Veq;
        weq_all(j,r) = weq;

        if all(real(lambda)<0)

            stability_all(j,r) = 1;

        elseif any(real(lambda)>0) && ...
               any(real(lambda)<0)

            stability_all(j,r) = -1;

        elseif all(real(lambda)>0)

            stability_all(j,r) = 0;

        else

            stability_all(j,r) = NaN;

        end

    end

end

%% Plot equilibrium stability

figure('Name','Question 15 - Equilibria');

hold on;

I_plot = repmat(I_values(:),1,max_eq);

stable = stability_all == 1;
saddle = stability_all == -1;
unstable = stability_all == 0;

h1 = plot(I_plot(stable),Veq_all(stable), ...
    'o','MarkerSize',5,'LineWidth',1.5);

h2 = plot(I_plot(saddle),Veq_all(saddle), ...
    'x','MarkerSize',7,'LineWidth',1.5);

h3 = plot(I_plot(unstable),Veq_all(unstable), ...
    's','MarkerSize',5,'LineWidth',1.5);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('V_{eq} (mV)');

title('Question 15: Equilibrium Points and Stability');

grid on;

legend([h1 h2 h3], ...
       {'Stable','Saddle','Unstable'}, ...
       'Location','best');

hold off;

%% Plot w equilibrium

figure('Name','Question 15 - w Equilibrium');

hold on;

plot(I_values,weq_all(:,1),'LineWidth',2);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('w_{eq}');

title('Question 15: Equilibrium w vs Applied Current');

grid on;

hold off;


%% ========================================================================
% QUESTION 16
% Reduced Hodgkin-Huxley Model
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 16\n');
fprintf('============================================================\n');

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

alpha_n = @(V) alpha_n_safe(V);
beta_n = @(V) 0.125*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V);
beta_m = @(V) 4*exp(-(V+60)/18);

alpha_h = @(V) 0.07*exp(-(V+60)/20);
beta_h = @(V) 1./(1+exp(-(V+30)/10));

%% Steady-state m and h

m_inf = @(V) ...
    alpha_m(V)./(alpha_m(V)+beta_m(V));

h_inf = @(V) ...
    alpha_h(V)./(alpha_h(V)+beta_h(V));

%% Initial condition

V0 = -60;

n0 = alpha_n(V0)/(alpha_n(V0)+beta_n(V0));

Y0 = [V0;n0];

%% Reduced HH

MLE = @(t,Y) [ ...
    (Iext ...
    -GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*m_inf(Y(1))^3*h_inf(Y(1))*(Y(1)-ENa) ...
    -GL*(Y(1)-EL))/C;

    phi*( ...
    alpha_n(Y(1))*(1-Y(2)) ...
    -beta_n(Y(1))*Y(2))
];

%% Solve

tspan = [0 100];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

[t,Y] = ode15s(MLE,tspan,Y0,options);

V = Y(:,1);
n = Y(:,2);

%% Voltage

figure('Name','Question 16 - Voltage');

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');

title('Question 16: Reduced HH Model - Membrane Potential');

grid on;

%% n

figure('Name','Question 16 - n');

plot(t,n,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('n');

title('Question 16: Reduced HH Model - Potassium Activation');

grid on;


%% ========================================================================
% QUESTION 17
% Anode-Break Excitation
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 17\n');
fprintf('============================================================\n');

%% Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;

T = 6.3;
Q = 3;

phi = Q^((T-6.3)/10);

%% Resting state

Vrest = -60;
Irest = 0;

ninf = alpha_n_safe(Vrest) / ...
       (alpha_n_safe(Vrest)+beta_n(Vrest));

minf = alpha_m_safe(Vrest) / ...
       (alpha_m_safe(Vrest)+beta_m(Vrest));

hinf = alpha_h(Vrest) / ...
       (alpha_h(Vrest)+beta_h(Vrest));

%% Calculate EL

EL = Vrest - ...
     (-GK*ninf^4*(Vrest-EK) ...
      -GNa*minf^3*hinf*(Vrest-ENa) ...
      + Irest)/GL;

fprintf('\nResting-state values:\n');
fprintf('Vrest = %.6f mV\n',Vrest);
fprintf('ninf = %.6f\n',ninf);
fprintf('minf = %.6f\n',minf);
fprintf('hinf = %.6f\n',hinf);
fprintf('EL = %.6f mV\n',EL);

%% Anodal current

Iext = @(t) -3*(t<=20);

%% Initial condition

Y0 = [Vrest;ninf;minf;hinf];

%% HH equations

HH = @(t,Y) [ ...
    (Iext(t) ...
    -GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    -GL*(Y(1)-EL))/C;

    phi*( ...
    alpha_n_safe(Y(1))*(1-Y(2)) ...
    -beta_n(Y(1))*Y(2));

    phi*( ...
    alpha_m_safe(Y(1))*(1-Y(3)) ...
    -beta_m(Y(1))*Y(3));

    phi*( ...
    alpha_h_safe(Y(1))*(1-Y(4)) ...
    -beta_h(Y(1))*Y(4))
];

%% Simulation

tspan = [0 1000];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

[t,Y] = ode15s(HH,tspan,Y0,options);

V = Y(:,1);
n = Y(:,2);
m = Y(:,3);
h = Y(:,4);

%% Voltage

figure('Name','Question 17 - Voltage');

plot(t,V,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');

title('Question 17: Anode-Break Excitation');

grid on;

%% Gating variables

figure('Name','Question 17 - Gating');

plot(t,n,'LineWidth',1.5);
hold on;

plot(t,m,'LineWidth',1.5);
plot(t,h,'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('Gating variable');

title('Question 17: HH Gating Variables');

legend('n','m','h');

grid on;

hold off;


%% ========================================================================
% QUESTION 18
% V-m Phase Plane with n,h Fixed
% ========================================================================

fprintf('\n============================================================\n');
fprintf('QUESTION 18\n');
fprintf('============================================================\n');

%% Parameters

C = 1;

GK  = 36;
GNa = 120;
GL  = 0.3;

EK  = -72;
ENa = 55;
EL  = -49.4;

phi = 1;

%% Rate functions

alpha_n = @(V) alpha_n_safe(V);
beta_n  = @(V) 0.125*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V);
beta_m  = @(V) 4*exp(-(V+60)/18);

alpha_h = @(V) 0.07*exp(-(V+60)/20);
beta_h  = @(V) 1./(exp(-(V+30)/10)+1);

%% Steady-state gating variables

ninf = @(V) ...
    alpha_n(V)./(alpha_n(V)+beta_n(V));

minf = @(V) ...
    alpha_m(V)./(alpha_m(V)+beta_m(V));

hinf = @(V) ...
    alpha_h(V)./(alpha_h(V)+beta_h(V));

%% Voltage and m ranges

V = linspace(-100,60,500);
m = linspace(0,1,500);

[VV,MM] = meshgrid(V,m);

%% CASE 1
% n and h fixed at resting values

Vrest = -60;

n_fixed1 = ninf(Vrest);
h_fixed1 = hinf(Vrest);

Vdot1 = ( ...
    -GK*n_fixed1^4.*(VV-EK) ...
    -GNa*MM.^3*h_fixed1.*(VV-ENa) ...
    -GL*(VV-EL))/C;

Mdot1 = ...
    alpha_m(VV).*(1-MM) ...
    -beta_m(VV).*MM;

figure('Name','Question 18 - Case 1');

hold on;

contour(VV,MM,Vdot1,[0 0], ...
    'LineWidth',2);

contour(VV,MM,Mdot1,[0 0], ...
    'LineWidth',2);

xlabel('V (mV)');
ylabel('m');

title('Question 18: V-m Phase Plane, n,h at Rest');

legend('V-nullcline','m-nullcline');

grid on;

xlim([-100 60]);
ylim([0 1]);

%% Find equilibria

F1 = @(Veq) ...
    -GK*n_fixed1^4*(Veq-EK) ...
    -GNa*minf(Veq)^3*h_fixed1*(Veq-ENa) ...
    -GL*(Veq-EL);

roots1 = [];

Vscan = linspace(-100,60,3000);

Fscan = arrayfun(F1,Vscan);

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1)<0

        r = fzero(F1,[Vscan(i),Vscan(i+1)]);

        if isempty(roots1) || ...
                all(abs(roots1-r)>1e-5)

            roots1(end+1) = r;

        end

    end

end

mroots1 = minf(roots1);

plot(roots1,mroots1,'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

fprintf('\nCASE 1: n,h at rest\n');

fprintf('n_fixed = %.6f\n',n_fixed1);
fprintf('h_fixed = %.6f\n',h_fixed1);

fprintf('Equilibrium points:\n');

for i = 1:length(roots1)

    fprintf('V = %.6f mV, m = %.6f\n', ...
        roots1(i),mroots1(i));

end

hold off;

%% CASE 2
% n and h at end of anodal stimulus

Iext = @(t) -3*(t<=20);

Y0 = [ ...
    Vrest;
    ninf(Vrest);
    minf(Vrest);
    hinf(Vrest)
];

HH = @(t,Y) [ ...
    (Iext(t) ...
    -GK*Y(2)^4*(Y(1)-EK) ...
    -GNa*Y(3)^3*Y(4)*(Y(1)-ENa) ...
    -GL*(Y(1)-EL))/C;

    alpha_n(Y(1))*(1-Y(2)) ...
    -beta_n(Y(1))*Y(2);

    alpha_m(Y(1))*(1-Y(3)) ...
    -beta_m(Y(1))*Y(3);

    alpha_h(Y(1))*(1-Y(4)) ...
    -beta_h(Y(1))*Y(4)
];

options = odeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-10);

[t,Y] = ode15s(HH,[0 20],Y0,options);

n_fixed2 = Y(end,2);
h_fixed2 = Y(end,4);

fprintf('\nCASE 2: n,h at end of anodal stimulus\n');

fprintf('n_fixed = %.6f\n',n_fixed2);
fprintf('h_fixed = %.6f\n',h_fixed2);

%% Case 2 phase plane

Vdot2 = ( ...
    -GK*n_fixed2^4.*(VV-EK) ...
    -GNa*MM.^3*h_fixed2.*(VV-ENa) ...
    -GL*(VV-EL))/C;

Mdot2 = ...
    alpha_m(VV).*(1-MM) ...
    -beta_m(VV).*MM;

figure('Name','Question 18 - Case 2');

hold on;

contour(VV,MM,Vdot2,[0 0], ...
    'LineWidth',2);

contour(VV,MM,Mdot2,[0 0], ...
    'LineWidth',2);

xlabel('V (mV)');
ylabel('m');

title('Question 18: V-m Phase Plane, n,h After Anodal Stimulus');

legend('V-nullcline','m-nullcline');

grid on;

xlim([-100 60]);
ylim([0 1]);

%% Case 2 equilibria

F2 = @(Veq) ...
    -GK*n_fixed2^4*(Veq-EK) ...
    -GNa*minf(Veq)^3*h_fixed2*(Veq-ENa) ...
    -GL*(Veq-EL);

roots2 = [];

Fscan = arrayfun(F2,Vscan);

for i = 1:length(Vscan)-1

    if Fscan(i)*Fscan(i+1)<0

        r = fzero(F2,[Vscan(i),Vscan(i+1)]);

        if isempty(roots2) || ...
                all(abs(roots2-r)>1e-5)

            roots2(end+1) = r;

        end

    end

end

mroots2 = minf(roots2);

plot(roots2,mroots2,'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

fprintf('\nCASE 2 equilibrium points:\n');

for i = 1:length(roots2)

    fprintf('V = %.6f mV, m = %.6f\n', ...
        roots2(i),mroots2(i));

end

hold off;


%% ========================================================================
% LOCAL FUNCTIONS
% IMPORTANT:
% ALL LOCAL FUNCTIONS MUST BE AT THE VERY END OF THE SCRIPT.
% ========================================================================

%% Safe alpha_n

function a = alpha_n_safe(V)

    a = -0.01.*(V+50) ./ ...
        (exp(-(V+50)/10)-1);

    idx = abs(V+50)<1e-8;

    a(idx) = 0.1;

end


%% Safe alpha_m

function a = alpha_m_safe(V)

    a = -0.1.*(V+35) ./ ...
        (exp(-(V+35)/10)-1);

    idx = abs(V+35)<1e-8;

    a(idx) = 1;

end


%% beta_n

function b = beta_n(V)

    b = 0.125.*exp(-(V+60)/80);

end


%% beta_m

function b = beta_m(V)

    b = 4.*exp(-(V+60)/18);

end


%% alpha_h

function a = alpha_h_safe(V)

    a = 0.07.*exp(-(V+60)/20);

end


%% beta_h

function b = beta_h_safe(V)

    b = 1./(exp(-(V+30)/10)+1);

end

