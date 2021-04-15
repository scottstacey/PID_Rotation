function RhaS
clear; clc; clf;
% This is a model of a PI controller using RhaS as the ouput
% Integral control provided via sRNA.
%% Initialising Global Parameters 

global k_1       % Maximal transcription rate of Z1
global U         % Concentration of Inducer U in nM
global K_U       % Concentraton of U at which Z1 transscription is at half maximal rate nM
global gamma_m   % Degrgadation and dilution rate of mRNA 
global theta     % Annihilation/binding rate of sRNA and mRNA 
global k_2       % Translation rate of RhaS from Z1 mRNA
global delta_1   % Degradation and diltuion rate of RhaS
global alpha_1   % Dimerisation rate of RhaS
global alpha_2   % Dissociation rate of RhaS dimer
global delta_2   % Degradation and dilution rate of RhaS Dimer
global k_3       % Maximal rate of transcription of Z2 sRNA
global K_X       % Concentration of X (RhaS Dimer) at which Z2 is transcribed at half maximal rate nM
global gamma_s   % Degradation and dilution rate of sRNA 
global setpoint  % The setpoint concentration of the output X under idealised integral control.

%% Parameter Values

k_1      = 0.1 * 60 * 60;
U        = 250000;
K_U      = 178000; 
theta    = 0.025 * 60 * 60;
k_2      = 0.06 * 60 * 60;
delta_1  = 0.00039 * 60 * 60;
alpha_1  = 0.2 * 60 * 60;
alpha_2  = 1 * 60 * 60;
delta_2  = delta_1;
k_3      = 1.5 * 60 * 60; 
K_X      = 2600; 
gamma_s  = 0.0008 * 60 * 60 * 0;
gamma_m  = 0.0041 * 60 * 60 * 0;

setpoint = (k_1 * K_X * U) / ((k_3 * K_U) + (k_3 * U) - (k_1 * U)); 

%% State is [Z1 Z2 RhaS X]
s0       = [0 0 0 0]; % Initial values of the states in the ODE model 
%% Generate the simulation 
Tend     = 15;        
ODEFUN   = @RhaSddt;
[t, S]   = ode45(ODEFUN, [0,Tend], s0);

%% Generate Plot Figure
figure(1);
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'b--', t, S(:,2), 'r', t, S(:,3), 'b', t, S(:,4), 'k', t, (2 * S(:,4) + S(:,3)), 'g', 'LineWidth', 3)
yline(setpoint, 'k--', 'LineWidth', 3);
legend('Z1 (mRNA)', 'Z2 sRNA', 'RhaS', 'RhaS Dimer', 'GFP', 'Setpoint', 'Location', 'northeast')
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS Control Circuit');


end

%% Dynamics 
function dS = RhaSddt(t, S);

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

Z_1  = S(1);
Z_2  = S(2);
RhaS = S(3);
X    = S(4);

dZ_1dt  = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt  = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dRhaSdt = (k_2 * Z_1) - (delta_1 * RhaS) - (2 * alpha_1 * RhaS * RhaS) + (2 * alpha_2 * X);
dXdt    = (alpha_1 * RhaS * RhaS) - (alpha_2 * X) - (delta_2 * X);

dS = [dZ_1dt; dZ_2dt; dRhaSdt; dXdt];

end