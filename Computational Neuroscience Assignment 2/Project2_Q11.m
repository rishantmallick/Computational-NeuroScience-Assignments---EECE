clear;
clc;
close all;
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
