clear; clc; clf;
%% Parameter Values
tau_1    = 0.09;        % Leaky transcription rate of Z_1 
k_1      = 0.85;        % Maximal transcription rate of Z_1
theta    = 0.05;        % Annihilation rate of Z_1 and Z_2 
opsilion = 0.0008;      % sRNA degradation/dilution rate
k_2      = 0.000000224; % Z_1 sRNA and X mRNA binding rate 
B        = 2600;        % A binding rate for U mediated induction of Z_1 
U        = 0.45;        % Concentration of the inducer for Z_1 transcription 
tau_2    = 0.09;        % Leaky transcription of Z_2 
k_3      = 0.85;        % Maximal transcription rate of 

%% 