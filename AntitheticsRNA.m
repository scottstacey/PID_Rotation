function AntitheticsRNA
clear; clc; clf;
% This is a model of the antithetic controller under translational control
% via sRNAs.
%% Initialising Global Parameters 
global tau_1
global k_1
global theta
global delta_s
global k_2
global B
global U 
global tau_2 
global k_3
global K
global alpha 
global delta_m 
global delta_c
global delta_p 
global T_1
global T_0 

%% Parameter Values
tau_1    = 0.09;        % Leaky transcription rate of Z_1 
k_1      = 0.85;        % Maximal transcription rate of Z_1
theta    = 0.05;        % Annihilation rate of Z_1 and Z_2 
delta_s  = 0.0008;      % sRNA degradation/dilution rate
k_2      = 0.000000224; % Z_1 sRNA and X mRNA binding rate 
B        = 2600;        % A binding rate for U mediated induction of Z_1 
U        = 0.45;        % Concentration of the inducer for Z_1 transcription 
tau_2    = 0.09;        % Leaky transcription of Z_2 
k_3      = 0.85;        % Maximal transcription rate of Z_2
K        = 600;         % Dissociation constant for X Protein (this is a pretty random number)
alpha    = 0.85;        % Transcription rate of X mRNA 
delta_m  = 0.0041;      % degradation rate of X mRNA 
delta_c  = 0.0041;      % degradation rate of sRNA-mRNA complex 
delta_p  = 0.00039;     % degradation rate of the X protein
T_1      = 0.0243;      % Translation rate of X from X mRNA 
T_0      = 0;           % Leaky translation rate of X from the sRNA-mRNA complex - currently assuming no leak

%% State is [Z_1 Z_2 X_m X_P C]
s0       = [0 0 0 0 0]; % Initial values of the states in the ODE model 

%% Generate the simulation 
Tend     = 2500;        % End time value -- This is currently random 
ODEFUN   = @antitheticsrnaddt;
[t, S] = ode45(ODEFUN, [0,Tend], s0);

%% Generate Plot Figure
figure(1);
set(gca, 'fontsize', 14);
plot(t, S(:,4), 'k', t, S(:,1), 'r', t, S(:,2), 'g', t, S(:,3), 'b', t, S(:,5), 'y', 'LineWidth', 3);
legend('X Protein', 'Z1', 'Z2', 'X mRNA', 'X-mRNA-Z1 Complex', 'Location', 'northwest');
xlabel('Time (seconds)');

end

%% Dynamics 
function dS = antitheticsrnaddt(t, S);

global tau_1;
global k_1;
global theta;
global delta_s;
global k_2;
global B;
global U; 
global tau_2; 
global k_3;
global K;
global alpha; 
global delta_m; 
global delta_c;
global delta_p;
global T_1;
global T_0;

Z_1 = S(1);
Z_2 = S(2);
X_m = S(3);
X_P = S(4);
C   = S(5);

dZ_1dt = tau_1 + ((k_1 * U)/(B + U)) - (theta * Z_1 * Z_2) - (delta_s * Z_1) - (k_2 * Z_1 * X_m);
dZ_2dt = tau_2 + ((k_3 * K) / (K + X_P)) - (delta_s * Z_2) - (theta * Z_1 * Z_2);
dX_mdt = alpha - (delta_m * X_m) - (k_2 * Z_1 * X_m);
dX_Pdt = (T_1 * X_m) + (T_0 * C) - (delta_p * X_P);
dCdt   = (k_2 * Z_1 * X_m) - (delta_c * C);

dS = [dZ_1dt; dZ_2dt; dX_mdt; dX_Pdt; dCdt];

end 