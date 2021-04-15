function RhaSPI1
clear; clc; clf;
% This is a model of a PI controller using RhaS as the ouput. Which tracks
% changes in the input U.
% Integral control provided via sRNA annihilation.
%% Initialising Global Parameters 

global k_1       % Maximal transcription rate of Z1
global U1        % Initial concentration of Inducer U in nM
global U2        % Second concentration of inducer U in nM
global U3        % Third concentration of inducer U in nM
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
global setpoint_1% Steady state value of X for U1
global setpoint_2% Steady state value of X for U2
global setpoint_3% Steady state value of X for U3

%% Parameter Values

k_1      = 0.1 * 60 * 60;
%k_1      = 0.85 * 60 * 60; 
K_U      = 178000; 
gamma_m  = 0.0041 * 60 * 60 * 0;
theta    = 0.025 * 60 * 60 ;
%theta    = 0.000000224 * 60 * 60;
k_2      = 0.06 * 60 * 60;
%k_2      = 0.0243 * 60 * 60;
alpha_1  = .1 * 60 * 60;
alpha_2  = .2 * 60 * 60;
k_3      = 1.5 * 60 * 60; 
%k_3      = 1670 * 60 * 60; 
K_X      = 2600; 
gamma_s  = 0.0008*60*60*0;
U1       = 250000;
U2       = 2500000;
U3       = 25000;
delta_1  = 0.039 * 60 * 60;
delta_2  = 0.039 * 60 * 60;
setpoint_1 = ((k_1 * K_X * U1) / ((k_3 * K_U) + (k_3 * U1) - (k_1 * U1)));
setpoint_2 = ((k_1 * K_X * U2) / ((k_3 * K_U) + (k_3 * U2) - (k_1 * U2)));
setpoint_3 = ((k_1 * K_X * U3) / ((k_3 * K_U) + (k_3 * U3) - (k_1 * U3)));

%% State is [Z1 Z2 RhaS X U]
s0       = [0 0 0 0 U1]; % Initial values of the states in the ODE model 
%% Generate the simulation 
T    = 15;        
ODEFUN   = @RhaSPIddt1;
[t, S]   = ode45(ODEFUN, [0,T], s0);

S_final = S(end,:);
T1        = 25;
s1        = [S_final(1) S_final(2) S_final(3) S_final(4) U2];
[t1, S1]  = ode45(ODEFUN, [T, T1], s1);

S1_final  = S1(end,:);
T2        = 40;
s2        = [S1_final(1) S1_final(2) S1_final(3) S1_final(4) U3];
[t2, S2]  = ode45(ODEFUN, [T1, T2], s2);
%% Generate Plot Figure
figure(1);
hold on;
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'b--', t1, S1(:,1), 'b--', t2, S2(:,1), 'b--', t, S(:,2), 'r', t1, S1(:,2), 'r', ...
    t2, S2(:,2), 'r', t, S(:,3), 'b', t1, S1(:,3), 'b', t2, S2(:,3), 'b', t, S(:,4), 'k', ...
    t1, S1(:,4), 'k', t2, S2(:,4), 'k', t, (2 * S(:,4) + S(:,3)), 'g', t1, (2 * S1(:,4) + S1(:,3)), 'g', ...
    t2, (2 * S2(:,4) + S2(:,3)), 'g', 'LineWidth', 3)
% plot([0 T], [setpoint_1 setpoint_1], 'k--', 'LineWidth', 3)
% plot([T T], [setpoint_1 setpoint_2], 'k--', 'LineWidth', 3)
% plot([T T1], [setpoint_2 setpoint_2], 'k--', 'LineWidth', 3)
% plot([T1 T1], [setpoint_2 setpoint_3], 'k--', 'LineWidth', 3)
% plot([T1 T2], [setpoint_3 setpoint_3], 'k--', 'LineWidth', 3)

h = zeros(6, 1);
h(1) = plot(NaN,NaN,'b--', 'LineWidth', 3);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 3);
h(3) = plot(NaN,NaN,'b', 'LineWidth', 3);
h(4) = plot(NaN,NaN,'k', 'LineWidth', 3);
h(5) = plot(NaN,NaN,'g', 'LineWidth', 3);
h(6) = plot(NaN,NaN,'k--', 'LineWidth', 3);
legend(h, 'Z_1 (mRNA)', 'Z_2 (sRNA)', 'RhaS', 'RhaS Dimer (X)', 'GFP', 'Location', 'northeast', 'NumColumns',2, 'LineWidth', 3);

xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('RhaS PI Controller');
hold off;


end

%% Dynamics 
function dS = RhaSPIddt1(t, S);

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

Z_1       = S(1);
Z_2       = S(2);
RhaS      = S(3);
X         = S(4);
U         = S(5);

dZ_1dt  = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt  = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dRhaSdt = (k_2 * Z_1) - (delta_1 * RhaS) - (2 * alpha_1 * RhaS * RhaS) + (2 * alpha_2 * X);
dXdt    = (alpha_1 * RhaS * RhaS) - (alpha_2 * X) - (delta_2 * X);
dUdt    = 0;

dS = [dZ_1dt; dZ_2dt; dRhaSdt; dXdt; dUdt];

end
