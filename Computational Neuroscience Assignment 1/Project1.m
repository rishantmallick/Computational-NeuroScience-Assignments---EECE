clc;
close all;
clear;
function dxdt = vanderpol(t,x,mu)
dxdt = zeros(2,1);
dxdt(1) = mu*x(2);
dxdt(2) = mu*(1-x(1)^2)*x(2) - x(1)/mu;  
end



%Problem 2 

x0 = [1;0];
mu_values = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
for i = 1:length(mu_values)

    mu = mu_values(i);
 
    tspan = [0,200];

    fprintf('Current value of mu = %.3f\n', mu);

    [t,x] = ode45(@(t,x) vanderpol(t,x,mu), tspan, x0);
  

    figure;                      % Creates a new figure window
    plot(t, x(:,1), 'LineWidth', 2);
    ylim([-4 4])
    title(['Van der Pol Oscillator through ode45 (\mu = ', num2str(mu), ')']);

    xlim([0 200]); 

    xlabel('Time');
    ylabel('y(t)');
    grid on;
    
end

%Problem 3

x0 = [1;0];
mu_values = [0.1,1,100,1000];
ode45_time = zeros(length(mu_values),1);
ode15s_time = zeros(length(mu_values),1);
for i = 1:length(mu_values)

    mu = mu_values(i);
    if mu==0.1 || mu==1 
      tspan = [0,200];
    elseif mu==100
      tspan = [0 3000];
    else
       tspan = [0,40000]; 
    end
    fprintf('Current value of mu = %.1f\n', mu);

    % ode45
    tic
    [t,x] = ode45(@(t,x) vanderpol(t,x,mu), tspan, x0);
    ode45_time(i) = toc;

    figure;                      % Creates a new figure window
    plot(t, x(:,1), 'LineWidth', 2);
    ylim([-4 4])
    title(['Van der Pol Oscillator through ode45 (\mu = ', num2str(mu), ')']);
    fprintf('Current value of ode45 = %.1f\n', ode45_time(i));
    if mu==0.1 || mu==1 
      xlim([0 200]);
    elseif mu==100
        xlim([0 3000]);
    else
       xlim([0 40000]); 
    end
    xlabel('Time');
    ylabel('y(t)');
    grid on;
    
end
for i = 1:length(mu_values)

    mu = mu_values(i);
    if mu==0.1 || mu==1 
      tspan = [0,200];
    elseif mu==100
      tspan = [0 3000];
    else
       tspan = [0,40000]; 
    end
    fprintf('Current value of mu = %.1f\n', mu);

    tic
    [t,x] = ode15s(@(t,x) vanderpol(t,x,mu), tspan, x0);
    ode15s_time(i) = toc;
    figure;                      % Creates a new figure window
    plot(t, x(:,1), 'LineWidth', 2);
    ylim([-4 4])
    fprintf('Current value of ode15s = %.1f\n', ode15s_time(i));
    title(['Van der Pol Oscillator through ode15s (\mu = ', num2str(mu), ')']);
    if mu==0.1 || mu==1 
      xlim([0 200]);
    elseif mu==100
        xlim([0 3000]);
    else
       xlim([0 40000]); 
    end
    xlabel('Time');
    ylabel('y(t)');
    grid on;
    
end
results = table(mu_values', ode45_time, ode15s_time,'VariableNames', {'Mu','ODE45_Time','ODE15s_Time'});
writetable(results, 'comparison.csv');
%question4
mu = 1;
tspan = [0 1000];


IC = [ 1   0;
       3   0;
      -2   2];

figure;
hold on;
grid on;
box on;

colors = lines(size(IC,1));

for i = 1:size(IC,1)

    x0 = IC(i,:)';

    % Solve ODE
    if mu >= 100
        [t,x] = ode15s(@(t,x) vanderpol(t,x,mu), tspan, x0);
    else
        [t,x] = ode45(@(t,x) vanderpol(t,x,mu), tspan, x0);
    end

    % Phase trajectory
    plot(x(:,1),x(:,2),...
        'Color',colors(i,:),...
        'LineWidth',2,...
        'DisplayName',sprintf('IC = [%g,%g]',x0(1),x0(2)));

    % Initial point
    plot(x0(1),x0(2),...
        'o',...
        'Color',colors(i,:),...
        'MarkerFaceColor',colors(i,:),...
        'MarkerSize',8,...
        'HandleVisibility','off');

    % Final point (shows convergence)
    plot(x(end,1),x(end,2),...
        'ks',...
        'MarkerFaceColor','k',...
        'MarkerSize',7,...
        'HandleVisibility','off');
end

xlabel('y');
ylabel('(1/\mu) dy/dt');

title(sprintf('Phase Plane of Van der Pol Oscillator (\\mu = %g)',mu));

legend('Location','best');

axis equal;
xlim([-3 3]);
ylim([-3.5 3.5]);

hold off;

