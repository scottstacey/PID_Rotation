function RhaS1
clear; clc; clf;
% This is a more complex version of the RhaS PI controller where there is a
% second inducer U2 (L-rhamnose) controlling the activity of the RhaS itself. 
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
global K_U2      % Concentration of U2 inducer (L-rhamnose) at which half X is bound by U2 and therefore active
global U2        % Concentration of L-rhamnose
global setpoint  % The setpoint concentration of the output X under idealised integral control.

%% Parameter Values

%k_1      = 0.85 * 60 * 60; 
k_1      = 0.1 * 60 * 60;
U        = 1000000;
K_U      = 178000; 
gamma_m  = 0.0041*60*60*0;
%theta    = 0.000000224 * 60 * 60;
theta    = 0.05 * 60 ;
%k_2      = 0.0243 * 60 * 60;
k_2      = 0.06 * 60 * 60;
delta_1  = 0.00039*60*60;
alpha_1  = 1 * 60 * 60;
alpha_2  = 1 * 60 * 60;
delta_2  = 0.00039*60*60;
%k_3      = 1670 * 60 * 60; 
k_3      = 1.5 * 60 * 60; 
K_X      = 1036000;
K_U2     = 1400;
U2       = 1000000;
gamma_s  = 0.0008*60*60*0;
setpoint = -(K_X * k_1 * U * (K_U2 + U2))/(U2 * (-k_3 * K_U - k_3 * U + k_1 * U)); 

%% State is [Z1 Z2 RhaS X]
s0       = [0 0 0 0]; % Initial values of the states in the ODE model 
%% Generate the simulation 
Tend     = 15;        % End time value -- This is currently random 
ODEFUN   = @RhaS1ddt;
[t, S]   = ode45(ODEFUN, [0,Tend], s0);

%% Generate Plot Figure
figure(1);
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'b--', t, S(:,3), 'b', t, S(:,4), 'k', t, (2 * S(:,4) + S(:,3)), 'g', 'LineWidth', 3)
yline(setpoint, 'k--', 'LineWidth', 3);
legend('Z1 (mRNA)', 'RhaS', 'RhaS Dimer', 'GFP', 'Location', 'northeast')
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS PI Controller');

figure(2);
set(gca, 'fontsize', 12);
plot(t, S(:,2), 'r', 'LineWidth', 3)
legend('Z2 (sRNA)', 'Location', 'northeast')
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS PI Controller');

end

%% Dynamics 
function dS = RhaS1ddt(t, S);

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
global U2
global K_U2

Z_1  = S(1);
Z_2  = S(2);
RhaS = S(3);
X    = S(4);

dZ_1dt  = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt  = ((k_3 * ((X * U2)/(K_U2 + U2)))/(K_X + ((X * U2)/(K_U2 + U2)))) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dRhaSdt = (k_2 * Z_1) - (delta_1 * RhaS) - (2 * alpha_1 * RhaS * RhaS) + (2 * alpha_2 * X);
dXdt    = (alpha_1 * RhaS * RhaS) - (alpha_2 * X) - (delta_2 * X);

dS = [dZ_1dt; dZ_2dt; dRhaSdt; dXdt];

end
