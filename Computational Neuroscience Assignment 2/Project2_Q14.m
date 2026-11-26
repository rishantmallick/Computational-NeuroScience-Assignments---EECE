clear;
clc;
close all;
%%Question14
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

alpha_n = @(V) alpha_n_safe(V);

beta_n = @(V) ...
    0.125.*exp(-(V+60)/80);

alpha_m = @(V) alpha_m_safe(V);

beta_m = @(V) ...
    4.*exp(-(V+60)/18);

alpha_h = @(V) ...
    0.07.*exp(-(V+60)/20);

beta_h = @(V) ...
    1./(exp(-(V+30)/10)+1);

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

alpha_n_sym = -0.01*(Vs+50) / ...
              (exp(-(Vs+50)/10)-1);

beta_n_sym = 0.125*exp(-(Vs+60)/80);

alpha_m_sym = -0.1*(Vs+35) / ...
              (exp(-(Vs+35)/10)-1);

beta_m_sym = 4*exp(-(Vs+60)/18);

f2 = 1*(alpha_n_sym*(1-ns) ...
      - beta_n_sym*ns);

f3 = 1*(alpha_m_sym*(1-ms) ...
      - beta_m_sym*ms);

f4 = 1*(alpha_h(Vs)*(1-hs) ...
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

    (alpha_n(Y(1))*(1-Y(2)) ...
    - beta_n(Y(1))*Y(2));

    (alpha_m(Y(1))*(1-Y(3)) ...
    - beta_m(Y(1))*Y(3));

    (alpha_h(Y(1))*(1-Y(4)) ...
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

function a = alpha_n_safe(V)

    a = -0.01.*(V+50) ./ ...
        (exp(-(V+50)/10)-1);

    idx = abs(V+50) < 1e-8;

    a(idx) = 0.1;

end


%% =========================================================
% Safe alpha_m
% ==========================================================

function a = alpha_m_safe(V)

    a = -0.1.*(V+35) ./ ...
        (exp(-(V+35)/10)-1);

    idx = abs(V+35) < 1e-8;

    a(idx) = 1;

end