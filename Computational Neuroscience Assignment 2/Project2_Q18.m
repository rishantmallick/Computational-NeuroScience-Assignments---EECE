clear;
clc;
close all;
%%Question 18
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

alpha_n = @(V) alpha_n_safe(V);
beta_n  = @(V) 0.125*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V);
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

function a = alpha_n_safe(V)

    a = -0.01.*(V+50) ./ ...
        (exp(-(V+50)/10)-1);

    idx = abs(V+50) < 1e-8;

    a(idx) = 0.1;

end


function a = alpha_m_safe(V)

    a = -0.1.*(V+35) ./ ...
        (exp(-(V+35)/10)-1);

    idx = abs(V+35) < 1e-8;

    a(idx) =1;

end