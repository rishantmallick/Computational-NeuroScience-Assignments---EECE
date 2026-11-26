clear; 
clc; 
close all; 
fprintf('Question 9\n'); 
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
  %% MLE equations 
  Iext = I(i); 
  MLE = @(t,Y) [ ... 
    (-gca*minf(Y(1))*(Y(1)-vca) ... 
     -gk*Y(2)*(Y(1)-vk) ... 
     -gl*(Y(1)-vl) + Iext)/C; 
      
    phi*(winf(Y(1))-Y(2))/tauw(Y(1)) 
  ]; 

  F = @(V) -gca*minf(V).*(V-vca) ... 
         -gk*winf(V).*(V-vk) ... 
         -gl*(V-vl) + Iext; 
 
  Veq = fzero(F,-30); 
  weq = winf(Veq);  

  % Morris-Lecar equations 
  syms V w 
  f1 = (-gca*minf(V)*(V-vca)-gk*w*(V-vk)-gl*(V-vl) + Iext) / C; 
  f2 = phi*(winf(V)-w) / tauw(V); 

  % Obtain Jacobian 
  J = jacobian([f1; f2], [V w]); 

  % Evaluate Jacobian at equilibrium 
  J_eq = double(subs(J, [V w], [Veq weq])); 

  % Calculate eigenvalues 
  eigenvalues = eig(J_eq); 

  % Display results 
  fprintf('Jacobian at equilibrium: %g\n', Iext) 
  format long 
  disp(J_eq) 
  fprintf('Eigen Values at equilibrium: %g\n', Iext) 
  disp(eigenvalues) 

  tspan = [0 3000]; 
  Y0=[Veq,weq]; 
  options = odeset('RelTol',1e-8,'AbsTol',1e-10); 
  [t,Y] = ode45(MLE,tspan,Y0,options); 
  V1t = Y(:,1); 

  % 2. Find spike times 
  figure; 
  plot(t,V1t); 
  xlabel('Time (ms)'); 
  ylabel('Voltage (mV)'); 
 
  titleName = ['V vs t, Iext = ', num2str(Iext)]; 
  title(titleName); 
  grid on; 
 
  filename = ['V_vs_t_Iext_' num2str(Iext) '_Q9.png']; 
  exportgraphics(gcf, filename, 'Resolution', 300); 

  [pks,locs] = findpeaks(V1t,t,'MinPeakHeight',-50); 

  % 3. Calculate inter-spike intervals (ISI) 
  ISI = diff(locs); 

  % 4. Calculate firing rate 
  if ~isempty(ISI) 
    firing_rate = 1000 / mean(ISI); 
  else 
    firing_rate = 0; 
  end 

  firingrate(i) = firing_rate; 

  % Display results 
  disp('Spike times (ms):'); 
  disp(locs); 
  disp('ISI (ms):'); 
  disp(ISI); 
  fprintf('Firing rate = %.4f Hz\n',firing_rate); 
end 

figure; 
plot(I, firingrate, 'o-'); 
xlabel('I_{ext} (\muA/cm^2)'); 
ylabel('Firing Rate (Hz)'); 
title('Firing Rate vs Applied Current'); 
grid on; 

exportgraphics(gcf, 'Firing_Rate_vs_Applied_Current_Q9.png', 'Resolution', 300);
