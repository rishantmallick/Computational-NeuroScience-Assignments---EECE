clear;
clc;
close all;
%%Question13
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

ninf = alpha_n_safe(Vrest) / ...
       (alpha_n_safe(Vrest) + beta_n(Vrest));

minf = alpha_m_safe(Vrest) / ...
       (alpha_m_safe(Vrest) + beta_m(Vrest));

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
    alpha_n_safe(Y(1))*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2));


    % dm/dt

    phi*( ...
    alpha_m_safe(Y(1))*(1-Y(3)) ...
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

function a = alpha_n_safe(V)

    if abs(V + 50) < 1e-8

        % Limit as V -> -50
        a = 0.1;

    else

        a = -0.01*(V+50) / ...
            (exp(-(V+50)/10)-1);

    end

end


%% Safe alpha_m function

function a = alpha_m_safe(V)

    if abs(V + 35) < 1e-8

        % Limit as V -> -35
        a = 1;

    else

        a = -0.1*(V+35) / ...
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