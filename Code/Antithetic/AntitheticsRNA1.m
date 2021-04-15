function AntitheticsRNA1
clear; clc; clf;
% This is a model of the antithetic controller which tracks changes in
% input U
%% Initialising Global Parameters 

global U1;        % Inducer concentration in nM
global U2;        % Inducer concentration in nM
global U3;        % Inducer concentration in nM
global k_1;      % Maximal rate of of Z1 (mRNA) transcription 
global K_U;      % A binding constant for U mediated induction of Z1 (mRNA) transcription
global gamma_m;  % Degradation rate of mRNA (Z1)
global gamma_s   % Degradation rate of sRNA (Z2)
global theta;    % Annihilation rate of Z1 (mRNA) and Z2 (sRNA)
global k_3;      % Maximal rate of Z2 (sRNA) transcription
global K_X;      % A binding constant for X mediatted induction of Z2 (sRNA) transcription
global k_2;      % Translation rate of X protein
global delta;    % Deradation/dilution rate of X (protein)
global setpoint_1% Steady state value of X for U1
global setpoint_2% Steady state value of X for U2
global setpoint_3% Steady state value of X for U3

%% Parameter Values

k_1      = 1*60*60;
k_2      = 0.06*60*60;
k_3      = 1.5*60*60;
theta    = 0.025 * 60 * 60;        % Same value as Khammash paper
gamma_m  = 0.0041*60*60*0;   % Set to zero to give integral control. Value taken from Kelly et al NAR
gamma_s  = 0.0008*60*60*0;   % Set to zero to give integral control. Value taken from Kelly et al NAR
K_U      = 178000;           %
K_X      = 2600;             %
delta    = 0.00039*60*60;    % Value taken from Kelly et al NAR
U1       = 250000;
U2       = 2500000;
U3       = 25000;
setpoint_1 = ((k_1 * K_X * U1) / ((k_3 * K_U) + (k_3 * U1) - (k_1 * U1)));
setpoint_2 = ((k_1 * K_X * U2) / ((k_3 * K_U) + (k_3 * U2) - (k_1 * U2)));
setpoint_3 = ((k_1 * K_X * U3) / ((k_3 * K_U) + (k_3 * U3) - (k_1 * U3)));

%% State is [Z1 Z2 X U]
s0       = [0 0 0 U1]; % Initial values of the states in the ODE model 
%% Generate the simulation 
T    = 15;        
ODEFUN   = @antitheticsrnaddt1;
[t, S]   = ode45(ODEFUN, [0,T], s0);

S_final = S(end,:);
T1        = 30;
s1        = [S_final(1) S_final(2) S_final(3) U2];
[t1, S1]  = ode45(ODEFUN, [T, T1], s1);

S1_final  = S1(end,:);
T2        = 50;
s2        = [S1_final(1) S1_final(2) S1_final(3) U3];
[t2, S2]  = ode45(ODEFUN, [T1, T2], s2);
%% Generate Plot Figure
figure(1);
hold on;
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'g', t1, S1(:,1), 'g', t2, S2(:,1), 'g', t, S(:,2), 'r', t1, S1(:,2), 'r', t2, S2(:,2), 'r', t, S(:,3), 'b', t1, S1(:,3), 'b', t2, S2(:,3), 'b', 'LineWidth', 3)
plot([0 T], [setpoint_1 setpoint_1], 'k--', 'LineWidth', 3)
plot([T T], [setpoint_1 setpoint_2], 'k--', 'LineWidth', 3)
plot([T T1], [setpoint_2 setpoint_2], 'k--', 'LineWidth', 3)
plot([T1 T1], [setpoint_2 setpoint_3], 'k--', 'LineWidth', 3)
plot([T1 T2], [setpoint_3 setpoint_3], 'k--', 'LineWidth', 3)

h = zeros(4, 1);
h(1) = plot(NaN,NaN,'g', 'LineWidth', 3);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 3);
h(3) = plot(NaN,NaN,'b', 'LineWidth', 3);
h(4) = plot(NaN,NaN,'k--', 'LineWidth', 3);
legend(h, 'Z_1 (mRNA)', 'Z_2 (sRNA)', 'X', 'Setpoint', 'Location', 'northeast', 'LineWidth', 3);

xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('Generic Antithetic Controller (With Dilution/degradation of Z_1 and Z_2)');
hold off;


end

%% Dynamics 
function dS = antitheticsrnaddt1(t, S);

global k_1;      % Maximal rate of of Z1 (mRNA) transcription 
global K_U;      % A binding constant for U mediated induction of Z1 (mRNA) transcription
global gamma_m;  % Degradation rate of mRNA (Z1)
global gamma_s   % Degradation rate of sRNA (Z2)
global theta;    % Annihilation rate of Z1 (mRNA) and Z2 (sRNA)
global k_3;      % Maximal rate of Z2 (sRNA) transcription
global K_X;      % A binding constant for X mediatted induction of Z2 (sRNA) transcription
global k_2;      % Translation rate of X protein
global delta;    % Deradation/dilution rate of X (protein)

Z_1 = S(1);
Z_2 = S(2);
X = S(3);
U = S(4);


dZ_1dt = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dXdt   = (k_2 * Z_1) - (delta * X);
dUdt   = 0;

dS = [dZ_1dt; dZ_2dt; dXdt; dUdt];

end