function RhaSPID
clear; clc; clf;
% This is a model of a PID controller using RhaS as the ouput. A mutant
% RhaS that maintains dimerisation but is unable to bind DNA (RhaS*) acts
% as a sink for RhaS providing derivative control.
% Integral control provided via sRNA annihilation.
%% Initialising Global Parameters 

global k_1
global U
global K_U
global gamma_m
global theta
global k_2 
global delta_1
global alpha_1
global alpha_2 
global delta_2
global k_3 
global K_X
global gamma_s
global k_4
global alpha_3
global alpha_4
global delta_3
global delta_4
global setpoint

%% Parameter Values

k_1      = 0.1 * 60 * 60;
U        = 250000;
K_U      = 178000; 
gamma_m  = 0.0041*60*60*0;
theta    = 0.05 * 60 ;
k_2      = 0.06 * 60 * 60;
delta_1  = 0.00039*60*60;
alpha_1  = 1 * 60 * 60;
alpha_2  = 1 * 60 * 60;
delta_2  = 0.00039*60*60;
k_3      = 1.5 * 60 * 60; 
K_X      = 2600; 
gamma_s  = 0.0008*60*60*0;
k_4      = 3000;
alpha_3  = 1 * 60 * 60;
alpha_4  = 1 * 60 * 60;
delta_3  = 0.00039*60*60 + 0;
delta_4  = 0.00039*60*60 + 0;
setpoint = (k_1 * K_X * U) / ((k_3 * K_U) + (k_3 * U) - (k_1 * U)); 

%% State is [Z1 Z2 RhaS X RhaS* X*]
s0       = [0 0 0 0 0 0]; % Initial values of the states in the ODE model 
%% Generate the simulation 
Tend     = 15;        % End time value -- This is currently random 
ODEFUN   = @RhaSPIDddt;
[t, S]   = ode45(ODEFUN, [0,Tend], s0);

%% Generate Plot Figure
figure(1);
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'b--', t, S(:,2), 'r', t, S(:,3), 'b', t, S(:,4), 'k', t, S(:,5), 'b:', t, S(:,6), 'k:', t, (2 * S(:,4) + S(:,3) + S(:,6)), 'g', 'LineWidth', 3)
yline(setpoint, 'k--', 'LineWidth', 3);
legend('Z1 (mRNA)', 'Z2 (sRNA)', 'RhaS', 'RhaS Dimer (X)', 'RhaS*', 'RhaS-RhaS* Dimer (X*)', 'GFP', 'Setpoint', 'Location', 'northeast')
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS PID Controller');


end

%% Dynamics 
function dS = RhaSPIDddt(t, S);

global k_1
global U
global K_U
global gamma_m
global theta
global k_2 
global delta_1
global alpha_1
global alpha_2 
global delta_2
global k_3 
global gamma_s
global K_X
global k_4
global alpha_3
global alpha_4
global delta_3
global delta_4

Z_1       = S(1);
Z_2       = S(2);
RhaS      = S(3);
X         = S(4);
RhaS_star = S(5);
X_star    = S(6);

dZ_1dt       = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt       = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dRhaSdt      = (k_2 * Z_1) - (delta_1 * RhaS) - (2 * alpha_1 * RhaS * RhaS) + (2 * alpha_2 * X) - (alpha_3 * RhaS * RhaS_star) + (alpha_4 * X_star);
dXdt         = (alpha_1 * RhaS * RhaS) - (alpha_2 * X) - (delta_2 * X);
dRhaS_stardt = k_4 - (delta_3 * RhaS_star) - (alpha_3 * RhaS * RhaS_star) + (alpha_4 * X_star);
dX_stardt    = (alpha_3 * RhaS * RhaS_star) - (alpha_4 * X_star) - (delta_4 * X_star);

dS = [dZ_1dt; dZ_2dt; dRhaSdt; dXdt; dRhaS_stardt; dX_stardt];

end
