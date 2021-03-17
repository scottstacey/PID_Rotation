function AntitheticsRNADamped
clear; clc; clf;
% This is a model of the antithetic controller under translational control
% via sRNAs.
%% Initialising Global Parameters 
global U;        % Inducer concentration in nM
global k_1;      % Maximal rate of of Z1 (mRNA) transcription 
global K_U;      % A binding constant for U mediated induction of Z1 (mRNA) transcription
global gamma_m;  % Degradation rate of mRNA (Z1)
global gamma_s   % Degradation rate of sRNA (Z2)
global theta;    % Annihilation rate of Z1 (mRNA) and Z2 (sRNA)
global k_3;      % Maximal rate of Z2 (sRNA) transcription
global K_X;      % A binding constant for X mediatted induction of Z2 (sRNA) transcription
global k_2;      % Translation rate of X protein
global delta;    % Deradation/dilution rate of X (protein)
global delta_c;  % Degradation rate of C (mRNA-sRNA complex)
global tau_1;    % Basal transcription rate of Z1 (X mRNA)
global tau_2;    % Basal transcription rate of Z2 (sRNA)
global tau_3;    % Leaky transaltion rate from C (mRNA-sRNA complex)
global setpoint; % Setpoint of steady state value of X
global alpha_1
global alpha_2
global delta_star

%% Parameter Values
% Parameter values in NM, s^-1, or NMs^-1.Multiplied by 60^2 to put in
% hours.
U        = 250000; 
k_1      = 0.1*60*60;
k_2      = 0.06*60*60;
k_3      = 1.5*60*60;
theta    = 0.05 * 60;        % Same value as Khammash paper, not sure where to find mRNA sRNA binding rates
gamma_m  = 0.0041*60*60*0;   % Set to zero to give integral control. Value taken from Kelly et al NAR
gamma_s  = 0.0008*60*60*0;   % Set to zero to give integral control. Value taken from Kelly et al NAR
K_U      = 178000;           %
K_X      = 2600;             %
delta    = 0.00039*60*60;    % Value taken from Kelly et al NAR
delta_c  = 0.0041*60*60;     % Value taken from Kelly et al NAR
tau_1    = 0*60*60;          % Set to 0 for simplicity 
tau_2    = 0*60*60;          % Set to 0 for simplicity
tau_3    = 0*60*60;          % Set to 0 for simplicity
setpoint = (k_1 * K_X * U) / ((k_3 * K_U) + (k_3 * U) - (k_1 * U)); % Setpoint for when gammas are 0 
alpha_1  = 0.1*60*60;
alpha_2  = 0.1*60*60;
delta_star = 0.039*60*60;

%% State is [Z_1 Z_2 X C X_star]
s0       = [0 0 0 0 0]; % Initial values of the states in the ODE model 
%% Generate the simulation 
Tend     = 12;        % End time value -- This is currently random 
ODEFUN   = @antitheticsrnaddt;
[t, S] = ode45(ODEFUN, [0,Tend], s0);

%% Generate Plot Figure
figure(1);
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'g', t, S(:,2), 'r', t, S(:,3), 'b', t, S(:,5), 'y', 'LineWidth', 3)
yline(setpoint, 'k--', 'LineWidth', 3);
legend('Z1 (x mRNA)','Z2 (sRNA)', 'X', 'X*', 'setpoint', 'Location', 'northeast')
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('sRNA Implementation of Antithetic Integral Feedback Circuit with Damping');


end

%% Dynamics 
function dS = antitheticsrnaddt(t, S);

global U;        % Inducer concentration in nM
global k_1;      % Maximal rate of of Z1 (mRNA) transcription 
global K_U;      % A binding constant for U mediated induction of Z1 (mRNA) transcription
global gamma_m;  % Degradation rate of mRNA (Z1)
global gamma_s   % Degradation rate of sRNA (Z2)
global theta;    % Annihilation rate of Z1 (mRNA) and Z2 (sRNA)
global k_3;      % Maximal rate of Z2 (sRNA) transcription
global K_X;      % A binding constant for X mediatted induction of Z2 (sRNA) transcription
global k_2;      % Translation rate of X protein
global delta;    % Deradation/dilution rate of X (protein)
global delta_c;  % Degradation rate of C (mRNA-sRNA complex)
global tau_1;    % Basal transcription rate of Z1 (X mRNA)
global tau_2;    % Basal transcription rate of Z2 (sRNA)
global tau_3;    % Leaky transaltion rate from C (mRNA-sRNA complex)
global alpha_1
global alpha_2
global delta_star

Z_1 = S(1);
Z_2 = S(2);
X = S(3);
C = S(4);
X_star = S(5);


dZ_1dt = tau_1 + ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt = tau_2 + ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dXdt   = (tau_3 * C) + (k_2 * Z_1) - (delta * X) - (alpha_1 * X) + (alpha_2 * X_star);
dCdt   = (theta * Z_1 * Z_2) - (delta_c * C);
dX_stardt = (alpha_1 * X) - (alpha_2 * X_star) - (delta_star * X_star); 

dS = [dZ_1dt; dZ_2dt; dXdt; dCdt; dX_stardt];

end 