function RhaSPID
clear; clc; clf;
% This is a model of a PID controller using RhaS as the ouput. A mutant
% RhaS that maintains dimerisation but is unable to bind DNA (RhaS*) acts
% as a sink for RhaS providing derivative control.
% Integral control provided via sRNA annihilation.
%% Initialising Global Parameters 

global k_1       % Maximal transcription rate of Z1
global U         % Concentration of Inducer U in nM
global K_U       % Concentraton of U at which Z1 transscription is at half maximal rate nM
global gamma_m   % Degrgadation and dilution rate of mRNA 
global theta     % Annihilation/binding rate of sRNA and mRNA 
global k_2       % Translation rate of RhaS from Z1 mRNA 
global delta_1   % Diltuion rate of RhaS based on E. coli doubling time of 30 minutes
global alpha_1   % Dimerisation rate of RhaS
global alpha_2   % Dissociation rate of RhaS dimer 
global delta_2   % Dilution rate of RhaS Dimer
global k_3       % Maximal rate of transcription of Z2 sRNA 
global K_X       % Concentration of X (RhaS Dimer) at which Z2 is transcribed at half maximal rate nM
global gamma_s   % Degradation and dilution rate of sRNA 
global k_4       % Prodcution rate of the mutant RhaS (RhaS*)
global delta_3   % Degradation and Dilution rate of RhaS*
global setpoint  % The setpoint concentration of the output X under idealised integral control. (When gammas=0)

%% Parameter Values

k_1      = 0.1 * 60 * 60;
U        = 250000;
K_U      = 178000; 
gamma_m  = 0.0041*60*60*0;  % set to zero for integral control
theta    = 0.025 * 60 * 60 ;
k_2      = 0.06 * 60 * 60;
delta_1  = 0.0078*60*60;
alpha_1  = 0.15001 * 60 * 60;
alpha_2  = 0.7001 * 60 * 60;
delta_2  = 0.0039*60*60;
k_3      = 1.5 * 60 * 60; 
K_X      = 2600; 
gamma_s  = 0.0008*60*60*0; % set to zero for integral control
k_4      = 1 * 60 * 60*0;              % Set to zero for PI controller, otherwise non-zero
delta_3  = (delta_2+delta_1)/2;             
setpoint = (k_1 * K_X * U) / ((k_3 * K_U) + (k_3 * U) - (k_1 * U)); 

%% State is [Z1 Z2 RhaS X RhaS* X*]
s0       = [0 0 0 0 0 0]; % Initial values of the states in the ODE model 
%% Generate the simulation 
Tend     = 15;        
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
title('RhaS^* Controller');

figure(2);
set(gca, 'fontsize', 12);
plot(t, S(:,4), 'k', 'LineWidth', 3)
yline(setpoint, 'k--', 'LineWidth', 3);
legend('RhaS Dimer (X)', 'Setpoint', 'Location', 'northeast')
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS^* Controller');


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
global delta_3



Z_1       = S(1);
Z_2       = S(2);
RhaS      = S(3);
X         = S(4);
RhaS_star = S(5);
X_star    = S(6);

dZ_1dt       = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt       = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dRhaSdt      = (k_2 * Z_1) - (delta_1 * RhaS) - (2 * alpha_1 * RhaS * RhaS) + (2 * alpha_2 * X) - (alpha_1 * RhaS * RhaS_star) + (alpha_2 * X_star);
dXdt         = (alpha_1 * RhaS * RhaS) - (alpha_2 * X) - (delta_1 * X);
dRhaS_stardt = k_4 - (delta_2 * RhaS_star) - (alpha_1 * RhaS * RhaS_star) + (alpha_2 * X_star);
dX_stardt    = (alpha_1 * RhaS * RhaS_star) - (alpha_2 * X_star) - (delta_3 * X_star);

dS = [dZ_1dt; dZ_2dt; dRhaSdt; dXdt; dRhaS_stardt; dX_stardt];

end
