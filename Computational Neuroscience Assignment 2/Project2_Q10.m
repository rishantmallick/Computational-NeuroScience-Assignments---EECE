clear;
clc;
close all;
fprintf('Question 10\n');
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
%% Equilibrium condition
% At equilibrium:
% dV/dt = 0
% dw/dt = 0
%
% Therefore w = winf(V)
F = @(V) (-gca*minf(V).*(V-vca) ...
          -gk*winf(V).*(V-vk) ...
          -gl*(V-vl) + Iext);
%% Find the three equilibrium voltages
Vscan = linspace(-100,100,2000);
Fscan = F(Vscan);
rootsV = [];
for i = 1:length(Vscan)-1
    if Fscan(i)*Fscan(i+1)<0
        root = fzero(F,[Vscan(i),Vscan(i+1)]);
        if isempty(rootsV) || all(abs(rootsV-root) > 1e-5)
            rootsV(end+1) = root;
        end
    end
end
Veq = rootsV;
weq = winf(Veq);
disp('Equilibrium points:')
disp([Veq(:), weq(:)])
for k = 1:length(Veq)
   syms Vsym wsym
   msym = 0.5*(1 + tanh((Vsym-V1)/V2));
   winfsym = 0.5*(1 + tanh((Vsym-V3)/V4));
   tauwsym = 1/cosh((Vsym-V3)/(2*V4));
   f1 = (-gca*msym*(Vsym-vca) ...
      -gk*wsym*(Vsym-vk) ...
      -gl*(Vsym-vl) + Iext)/C;
   f2 = phi*(winfsym-wsym)/tauwsym;
   J = jacobian([f1; f2],[Vsym wsym]);
   J_eq = double(subs(J,{Vsym,wsym},{Veq(k),weq(k)}));
   lambda = eig(J_eq);
   fprintf('Equilibrium %d:\n',k);
   fprintf('Veq = %.4f, weq = %.4f\n',Veq(k),weq(k));
   disp('Jacobian:');
   disp(J_eq);
   disp('Eigenvalues:');
   disp(lambda);
   if all(real(lambda) < 0)
    disp('Stable')
   elseif all(real(lambda) > 0)
    disp('Unstable')
   elseif any(real(lambda) > 0) && any(real(lambda) < 0)
    disp('Saddle point')
   end
end

V = linspace(-100,100,2000);
% V-nullcline: dV/dt = 0
V_null = (-gca*minf(V).*(V-vca) ...
          -gl*(V-vl) + Iext) ./ (gk*(V-vk));

% w-nullcline: dw/dt = 0
w_null = winf(V);
for k = 1:length(Veq)

    J_eq = double(subs(J,{Vsym,wsym},{Veq(k),weq(k)}));
    lambda = eig(J_eq);

    if real(lambda(1))*real(lambda(2)) < 0
        saddle_index = k;
    end
end

Vs = Veq(saddle_index);
Ws = weq(saddle_index);
J_saddle = double(subs(J,{Vsym,wsym},{Vs,Ws}));
[VEC,D] = eig(J_saddle);
lambda = diag(D);
unstable_index = find(real(lambda) > 0);
stable_index   = find(real(lambda) < 0);

v_unstable = VEC(:,unstable_index);
v_stable   = VEC(:,stable_index);
odefun = @(t,y) [
    (-gca*minf(y(1))*(y(1)-vca) ...
     -gk*y(2)*(y(1)-vk) ...
     -gl*(y(1)-vl) + Iext)/C;

    phi*(winf(y(1))-y(2))/tauw(y(1))
];
epsilon = 1e-3;
y0_u1 = [Vs; Ws] + epsilon*v_unstable;
y0_u2 = [Vs; Ws] - epsilon*v_unstable;

[t_u1,Y_u1] = ode45(odefun,[0 50],y0_u1);
[t_u2,Y_u2] = ode45(odefun,[0 50],y0_u2);
y0_s1 = [Vs; Ws] + epsilon*v_stable;
y0_s2 = [Vs; Ws] - epsilon*v_stable;

[t_s1,Y_s1] = ode45(odefun,[0 -50],y0_s1);
[t_s2,Y_s2] = ode45(odefun,[0 -50],y0_s2);
%% Figure 1: Nullclines and equilibrium points

figure;
hold on;

plot(V,V_null,'LineWidth',2);
plot(V,w_null,'LineWidth',2);

plot(Veq,weq,'ko','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('w');
title('Morris-Lecar Nullclines and Equilibrium Points');

legend('V-nullcline', ...
       'w-nullcline', ...
       'Equilibrium points');

grid on;
xlim([-50 15]);
ylim([0 0.4]);

hold off;


%% Figure 2: Stable and unstable manifolds

figure;
hold on;

% Unstable manifold
hU = plot(Y_u1(:,1),Y_u1(:,2),'-','LineWidth',2.5);
plot(Y_u2(:,1),Y_u2(:,2),'-','LineWidth',2.5);

% Stable manifold
hS = plot(Y_s1(:,1),Y_s1(:,2),'--','LineWidth',2.5);
plot(Y_s2(:,1),Y_s2(:,2),'--','LineWidth',2.5);

% Saddle point
hE = plot(Vs,Ws,'ko','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('w');
title('Stable and Unstable Manifolds');

legend([hU hS hE], ...
       'Unstable manifold', ...
       'Stable manifold', ...
       'Saddle point');

grid on;
xlim([-50 15]);
ylim([0 0.4]);

hold off;