clear;
clc;
close all;
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
    -0.01*(V+50) ./ ...
    (exp(-(V+50)/10)-1);

% beta_n
beta_n = @(V) ...
    0.125*exp(-(V+60)/80);

% alpha_m
alpha_m = @(V) ...
    -0.1*(V+35) ./ ...
    (exp(-(V+35)/10)-1);

% beta_m
beta_m = @(V) ...
    4*exp(-(V+60)/18);

% alpha_h
alpha_h = @(V) ...
    0.07*exp(-(V+60)/20);

% beta_h
beta_h = @(V) ...
    1./(exp(-(V+30)/10)+1);


%% Initial voltage

V0 = -60;


%% Handle 0/0 cases

% alpha_n has 0/0 at V = -50 mV
%
% lim V->-50 alpha_n(V) = phi

if abs(V0 + 50) < 1e-10
    an0 = 0.1;
else
    an0 = alpha_n(V0);
end
% alpha_m has 0/0 at V = -35 mV
%
% lim V->-35 alpha_m(V) = phi

if abs(V0 + 35) < 1e-10
    am0 = 1;
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
    alpha_n_safe(Y(1))*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2) );


    % dm/dt

    phi*( ...
    alpha_m_safe(Y(1))*(1-Y(3)) ...
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

function a = alpha_n_safe(V)

    if abs(V + 50) < 1e-8

        % Limit at V = -50 mV
        a = 0.1;

    else

        a = -0.01*(V+50) / ...
            (exp(-(V+50)/10)-1);

    end

end


%% Safe alpha_m

function a = alpha_m_safe(V)

    if abs(V + 35) < 1e-8

        % Limit at V = -35 mV
        a = 1;

    else

        a = -0.1*(V+35) / ...
            (exp(-(V+35)/10)-1);

    end

end