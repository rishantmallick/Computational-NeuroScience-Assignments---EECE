clear;
clc;
close all;
%%Question15
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

plot(I_values,weq_all(:,1),'LineWidth',2);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('w_{eq}');

title('Equilibrium w versus Applied Current');

grid on;

hold off;