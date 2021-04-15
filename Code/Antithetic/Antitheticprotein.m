function Antitheticprotein
clear; clc; clf;
% This is a model of the antithetic controller under transcriptional control
% via protein.
%% Initialising Global Parameters 
global k_1     % Maximal transcription rate of Z1 mRNA
global k_2     % Maximal transcription rate of X mRNA
global k_3     % Maximal transcription rate of Z2 mRNA 
global tau_1   % Leaky transcription of Z1 mRNA
global tau_2   % Leaky transcription of X mRNA
global tau_3   % Leaky transcription of Z2 mRNA 
global delta_p % Degradation/dilution rate of X Protein 
global gamma   % Degradation/dilution rate of Z1 and Z2 Protein
global theta   % Annihilation rate of Z1 and Z2 Proteins
global K_XP    % Binding constant for X Protein to regulate Z2 transcription
global K_Z1P   % Binding constant for Z1 Protein to regulate X transcription 
global K_U     % Binding constant for U to regulate Z1 transcription 
global U       % Concentration of inducer 

%% Parameter Values
U       = 100000;
k_1     = 0.1*60*60;         
k_2     = 0.1*60*60;         
k_3     = 0.5*60*60;         
tau_1   = 0;                 % Leaky transcription set to 0 for simplicity
tau_2   = 0;                 % Leaky transcription set to 0 for simplicity
tau_3   = 0;                 % Leaky transcription set to 0 for simplicity
gamma   = 0.028 * 60 * 0;        % Set as dilution rate for E. coli. Set to zero here for integral control?
delta_p = (0.00039 * 60 * 60) + gamma; % Set as same rate in antitheticsRNA.m
theta   = 0.05 * 60;         % Same rate as in Khammash Paper and antitheticsRNA.m
K_XP    = 2600;              % Same as antitheticsRNA.m
K_Z1P   = 100;              % Same as antitheticsRNA.m
K_U     = 178000;            % Same as antitheticsRNA.m

%% State is [Z_1 Z_2 X]
s0       = [0 0 0]; % Initial values of the states in the ODE model 

%% Generate the simulation 
Tend     = 30;        
ODEFUN   = @antitheticproteinddt;
[t, S] = ode45(ODEFUN, [0,Tend], s0);

%% Generate Plot Figure
figure(1);
set(gca, 'fontsize', 14);
plot(t, S(:,1), 'g', t, S(:,2), 'r', t, S(:,3), 'b', 'LineWidth', 3);
legend('Z1', 'Z2', 'X', 'Location', 'northeast');
xlabel('Time (Hours)');

end

%% Dynamics 
function dS = antitheticproteinddt(t, S);

global k_1     % Maximal transcription rate of Z1 mRNA
global k_2     % Maximal transcription rate of X mRNA
global k_3     % Maximal transcription rate of Z2 mRNA 
global tau_1   % Leaky transcription of Z1 mRNA
global tau_2   % Leaky transcription of X mRNA
global tau_3   % Leaky transcription of Z2 mRNA 
global delta_p % Degradation/dilution rate of X Protein 
global gamma   % Degradation/dilution rate of Z1 and Z2 Protein
global theta   % Annihilation rate of Z1 and Z2 Proteins
global K_XP    % Binding constant for X Protein to regulate Z2 transcription
global K_Z1P   % Binding constant for Z1 Protein to regulate X transcription 
global K_U     % Binding constant for U to regulate Z1 transcription 
global U       % Concentration of inducer 

Z_1 = S(1);
Z_2 = S(2);
X  = S(3);


dZ_1dt = tau_1 + ((k_1 * U)/(K_U + U)) - (gamma * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt = tau_3 + ((k_3 * X)/(K_XP + X)) - (gamma * Z_2) - (theta * Z_1 * Z_2);
dX_dt  = tau_2 + ((k_2 * Z_1)/(K_Z1P + Z_1)) - (delta_p * X);



dS = [dZ_1dt; dZ_2dt; dX_dt];

end 