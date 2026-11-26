clear;
clc;
close all;
%%Question16
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