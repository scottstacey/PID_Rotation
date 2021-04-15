function AntitheticsRNA2
clear; clc; clf;
% This is a model of the antithetic controller which tracks changes in
% response to a disturbance in the process under control -- i.e. translation by 
% varying the k_2 paramater
%% Initialising Global Parameters 

global U;        % Inducer concentration in nM
global k_1;      % Maximal rate of of Z1 (mRNA) transcription 
global K_U;      % A binding constant for U mediated induction of Z1 (mRNA) transcription
global gamma_m;  % Degradation rate of mRNA (Z1)
global gamma_s   % Degradation rate of sRNA (Z2)
global theta;    % Annihilation rate of Z1 (mRNA) and Z2 (sRNA)
global k_3;      % Maximal rate of Z2 (sRNA) transcription
global K_X;      % A binding constant for X mediatted induction of Z2 (sRNA) transcription
global k_2_1;    % Translation rate of X protein
global k_2_2;
global k_2_3;
global delta;    % Deradation/dilution rate of X (protein)
global setpoint% Steady state value of X for U1


%% Parameter Values

k_1      = 0.1*60*60;
k_2_1    = 0.06*60*60;
k_2_2    = 0.07*60*60;
k_2_3    = 0.05*60*60;
k_3      = 1.5*60*60;
theta    = 0.025 * 60 * 60;        % Same value as Khammash paper
gamma_m  = 0.0041*60*60*0;   % Set to zero to give integral control. Value taken from Kelly et al NAR
gamma_s  = 0.0008*60*60*0;   % Set to zero to give integral control. Value taken from Kelly et al NAR
K_U      = 178000;           %
K_X      = 2600;             %
delta    = 0.00039*60*60*100;    % Value taken from Kelly et al NAR
U        = 250000;
setpoint = ((k_1 * K_X * U) / ((k_3 * K_U) + (k_3 * U) - (k_1 * U)));

%% State is [Z1 Z2 X k_2]
s0       = [0 0 0 k_2_1]; % Initial values of the states in the ODE model 
%% Generate the simulation 
T    = 15;        
ODEFUN   = @antitheticsrnaddt2;
[t, S]   = ode45(ODEFUN, [0,T], s0);

S_final = S(end,:);
T1        = 30;
s1        = [S_final(1) S_final(2) S_final(3) k_2_2];
[t1, S1]  = ode45(ODEFUN, [T, T1], s1);

S1_final  = S1(end,:);
T2        = 50;
s2        = [S1_final(1) S1_final(2) S1_final(3) k_2_3];
[t2, S2]  = ode45(ODEFUN, [T1, T2], s2);
%% Generate Plot Figure
figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(figure(1),'defaultAxesColorOrder',[left_color; right_color]);
hold on;
set(gca, 'fontsize', 12);
plot(t, S(:,1), 'g', t1, S1(:,1), 'g', t2, S2(:,1), 'g', t, S(:,2), 'r', t1, S1(:,2), 'r', t2, S2(:,2), 'r', t, S(:,3), 'b', t1, S1(:,3), 'b', t2, S2(:,3), 'b', 'LineWidth', 3)
yline(setpoint, 'k--', 'LineWidth', 3);
plot(t, S(:,4), 'k', t1, S1(:,4), 'k', t2, S2(:,4), 'k', 'LineWidth', 3);
plot([T T], [k_2_1 k_2_2], 'k', 'LineWidth', 3);
plot([T1 T1], [k_2_2 k_2_3], 'k', 'LineWidth', 3);
ylabel('Concentration (nM)');
h = zeros(5, 1);
h(1) = plot(NaN,NaN,'g', 'LineWidth', 3);
h(2) = plot(NaN,NaN,'r', 'LineWidth', 3);
h(3) = plot(NaN,NaN,'b', 'LineWidth', 3);
h(4) = plot(NaN,NaN,'k--', 'LineWidth', 3);
h(5) = plot(NaN,NaN, 'k', 'LineWidth', 3);
legend(h, 'Z_1 (mRNA)', 'Z_2 (sRNA)', 'X', 'Setpoint', 'k_2 (Translation Rate)', 'Location', 'northeast', 'LineWidth', 3);


yyaxis right
ylim([0 350]);
ylabel('Translation Rate (k_2) (nM hr^{-1})');


xlabel('Time (Hours)');
title(["Response of RNA-based Antithetic Feedback Motif to Disturbances in k_2", "(With Dilution/Degradation of Z_1 or Z_2)"]);
hold off;


end

%% Dynamics 
function dS = antitheticsrnaddt2(t, S);

global k_1;      % Maximal rate of of Z1 (mRNA) transcription 
global U;
global K_U;      % A binding constant for U mediated induction of Z1 (mRNA) transcription
global gamma_m;  % Degradation rate of mRNA (Z1)
global gamma_s   % Degradation rate of sRNA (Z2)
global theta;    % Annihilation rate of Z1 (mRNA) and Z2 (sRNA)
global k_3;      % Maximal rate of Z2 (sRNA) transcription
global K_X;      % A binding constant for X mediatted induction of Z2 (sRNA) transcription
global delta;    % Deradation/dilution rate of X (protein)

Z_1 = S(1);
Z_2 = S(2);
X = S(3);
k_2 = S(4);


dZ_1dt = ((k_1 * U)/(K_U + U)) - (gamma_m * Z_1) - (theta * Z_1 * Z_2);
dZ_2dt = ((k_3 * X)/(K_X + X)) - (gamma_s * Z_2) - (theta * Z_1 * Z_2);
dXdt   = (k_2 * Z_1) - (delta * X);
dk_2dt = 0;

dS = [dZ_1dt; dZ_2dt; dXdt; dk_2dt];

end