function RhaSPIDS1
clear; clc; clf;
% This is a model of a PID controller using RhaS as the ouput. A mutant
% RhaS that maintains dimerisation but is unable to bind DNA (RhaS*) acts
% as a sink for RhaS providing derivative control.
% Integral control provided via sRNA annihilation.
%% Initialising Global Parameters 

global k_1       % Maximal transcription rate of Z1
global U1        % Initial concentration of Inducer U in nM
global U2        % Second concentration of inducer U in nM
global U3        % Third concentration of inducer U in nM
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
global alpha_3   % Dimmerisation rate between RhaS and RhaS*
global alpha_4   % Dissociation rate of X* (RhaS-RhaS*)
global deg       % Degradation rate of RhaS due to a degradation tag
global deg_star  % Degradation rate of RhaS* due to a degradation tag
global delta_3   % Degradation and Dilution rate of RhaS*
global delta_4   % Degradation and Dilution rate of X* 
global setpoint_1% Steady state value of X for U1
global setpoint_2% Steady state value of X for U2
global setpoint_3% Steady state value of X for U3

%% Parameter Values

k_1      = 0.1 * 60 * 60;
K_U      = 178000; 
gamma_m  = 0.0041*60*60*0;  % set to zero for integral control
theta    = 0.05 * 60 ;
k_2      = 0.06 * 60 * 60;
dilution = 0.00039*60*60;
deg      = 0.01*60*60;             % Degradation rate of RhaS due to a degradation tag
deg_star = 0.039*60*60;             % Degradation rate of RhaS* due to a degradation tag
delta_1  = dilution + deg;
alpha_1  = 1 * 60 * 60;
alpha_2  = 1 * 60 * 60;
delta_2  = delta_1;
k_3      = 1.5 * 60 * 60; 
K_X      = 2600; 
gamma_s  = 0.0008*60*60*0; % set to zero for integral control
k_4      = 15 * 60 * 60;              % Set to zero for PI controller, otherwise non-zero
alpha_3  = 1 * 60 * 60;
alpha_4  = 1 * 60 * 60;
delta_3  = dilution + deg_star; % Dilution Rate + Degradation Rate 
delta_4  = delta_3;             % Dilution Rate + Degrdation Rate 
U1       = 250000;
U2       = 2500000;
U3       = 25000;
setpoint_1 = ((k_1 * K_X * U1) / ((k_3 * K_U) + (k_3 * U1) - (k_1 * U1)));
setpoint_2 = ((k_1 * K_X * U2) / ((k_3 * K_U) + (k_3 * U2) - (k_1 * U2)));
setpoint_3 = ((k_1 * K_X * U3) / ((k_3 * K_U) + (k_3 * U3) - (k_1 * U3)));

%% State is [Z1 Z2 RhaS X RhaS* X*]
s0       = [0 0 0 0 0 0 U1]; % Initial values of the states in the ODE model 
%% Generate the simulation 
T    = 15;        
ODEFUN   = @RhaSPIDddt1;
[t, S]   = ode45(ODEFUN, [0,T], s0);

S_final = S(end,:);
T1        = 20;
s1        = [S_final(1) S_final(2) S_final(3) S_final(4) S_final(5) S_final(6) U2];
[t1, S1]  = ode45(ODEFUN, [T, T1], s1);

S1_final  = S1(end,:);
T2        = 25;
s2        = [S1_final(1) S1_final(2) S1_final(3) S1_final(4) S1_final(5) S1_final(6) U3];
[t2, S2]  = ode45(ODEFUN, [T1, T2], s2);
%% Generate Plot Figure
figure(1);
hold on;
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'b--', t1, S1(:,1), 'b--', t2, S2(:,1), 'b--', t, S(:,2), 'r', t1, S1(:,2), 'r', t2, S2(:,2), 'r', t, S(:,3), 'b', t1, S1(:,3), 'b', t2, S2(:,3), 'b', t, S(:,4), 'k', t1, S1(:,4), 'k', t2, S2(:,4), 'k', t, S(:,5), 'b:', t1, S1(:,5), 'b:', t2, S2(:,5), 'b:', t, S(:,6), 'k:', t1, S1(:,6), 'k:', t2, S2(:,6), 'k:', t, (2 * S(:,4) + S(:,3) + S(:,6)), 'g', t1, (2 * S1(:,4) + S1(:,3) + S1(:,6)), 'g', t2, (2 * S2(:,4) + S2(:,3) + S2(:,6)), 'g', 'LineWidth', 3)
plot([0 T], [setpoint_1 setpoint_1], 'k--', 'LineWidth', 3)
plot([T T], [setpoint_1 setpoint_2], 'k--', 'LineWidth', 3)
plot([T T1], [setpoint_2 setpoint_2], 'k--', 'LineWidth', 3)
plot([T1 T1], [setpoint_2 setpoint_3], 'k--', 'LineWidth', 3)
plot([T1 T2], [setpoint_3 setpoint_3], 'k--', 'LineWidth', 3)

h = zeros(7, 1);
h(1) = plot(NaN,NaN,'b--', 'LineWidth', 3);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 3);
h(3) = plot(NaN,NaN,'b', 'LineWidth', 3);
h(4) = plot(NaN,NaN,'k', 'LineWidth', 3);
h(5) = plot(NaN,NaN,'b:', 'LineWidth', 3);
h(6) = plot(NaN,NaN,'k:', 'LineWidth', 3);
h(7) = plot(NaN,NaN,'g', 'LineWidth', 3);
h(8) = plot(NaN,NaN,'k--', 'LineWidth', 3);
legend(h, 'Z1 (mRNA)', 'Z2 (sRNA)', 'RhaS', 'RhaS Dimer (X)', 'RhaS*', 'RhaS-RhaS* Dimer (X*)', 'GFP', 'Setpoint', 'Location', 'northwest', 'NumColumns',2, 'LineWidth', 3);

xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS PID Controller');
hold off;


end

%% Dynamics 
function dS = RhaSPIDddt1(t, S);

global k_1
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
global deg_star
global delta_3
global delta_4


Z_1       = S(1);
Z_2       = S(2);
RhaS      = S(3);
X         = S(4);
RhaS_star = S(5);
X_star    = S(6);
U         = S(7);

dZ_1dt       = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt       = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dRhaSdt      = (k_2 * Z_1) - (delta_1 * RhaS) - (2 * alpha_1 * RhaS * RhaS) + (2 * alpha_2 * X) - (alpha_3 * RhaS * RhaS_star) + (alpha_4 * X_star) + (deg_star * X_star);
dXdt         = (alpha_1 * RhaS * RhaS) - (alpha_2 * X) - (delta_2 * X);
dRhaS_stardt = k_4 - (delta_3 * RhaS_star) - (alpha_3 * RhaS * RhaS_star) + (alpha_4 * X_star);
dX_stardt    = (alpha_3 * RhaS * RhaS_star) - (alpha_4 * X_star) - (delta_4 * X_star);
dUdt         = 0;

dS = [dZ_1dt; dZ_2dt; dRhaSdt; dXdt; dRhaS_stardt; dX_stardt; dUdt];

end
