function AntitheticDamped
clear; clc; clf;
% This is a model of the antithetic controller under transcriptional control
% via protein.
%% Initialising Global Parameters 
global setpoint    % Setpoint of steady state value of X 
global mu          % Transcription rate of Z1 (actuator)
global eta         % Annihiliation Rate of Z1 and Z2
global theta       % Transcription rate constant of Z2 (Sensor)
global k           % Transcription rate constant of X
global gamma       % Degradation rate of X 
global delta       % Dilution rate
global gamma_star  % Degradation rate of X*
global alpha_1     % Production rate of X* from X
global alpha_2     % Production rate of X from X*

%% Parameter Values
% Parameters taken from the Khammash paper in units of nM and min^-1,
% multipled by 60 to give hr^-1. 
mu         = 0.028*60;         % Search space between 10^0 and 10^4 nM min^-1
eta        = 0.05*60;          % Conservative lower bound where higher eta values improve performance
theta      = 0.028*60;         % 0.028  min^-1 used as a reference for search space between 10^-2 and 10^2
setpoint   = mu / theta;       % Setpoint when delta = 0
k          = 0.028*60;         % 0.028  min^-1 used as a reference for search space between 10^-2 and 10^2
gamma      = 0.014*60;         % Native degradation of output set at 0.5 * delta
delta      = 0.028*60*0;       % Set to zero for the case of the idealised Antithetic integral feedback circuit 
gamma_star = 0.084*60;
alpha_1    = 0.1*60;
alpha_2    = 0.1*60;

%% State is [Z_1 Z_2 X X_star]
s0       = [0 0 0 0];  % Initial values of the states in the ODE model 
%% Generate the simulation 
Tend     = 30;       % End time value -- This is currently random 
ODEFUN   = @antitheticddt;
[t, S] = ode45(ODEFUN, [0,Tend], s0);


%% Generate Plot Figure
figure(1);
hold on;
set(gca, 'fontsize', 12);
plot(t, S(:,3), 'b', t, S(:,1),'g', t, S(:,2), 'r', t, S(:,4), 'y', 'LineWidth', 3);
yline(setpoint, 'k--', 'LineWidth', 3);
legend('X', 'Z1', 'Z2', 'X*', 'Set Point', 'Location', 'northeast');
xlabel('Time (Hours)');
ylabel('Concentration (nM)');
title('Idealised Antithetic Integral Feedbak Controller With Damping');
hold off;

end

%% Dynamics 
function dS = antitheticddt(t, S);

global mu     % Transcription rate of Z1 (actuator)
global eta    % Annihiliation Rate of Z1 and Z2
global theta  % Transcription rate constant of Z2 (Sensor)
global k      % Transcription rate constant of X
global gamma  % Degradation rate of X 
global delta  % Dilution rate
global gamma_star  % Degradation rate of X*
global alpha_1     % Production rate of X* from X
global alpha_2     % Production rate of X from X*

Z1     = S(1);
Z2     = S(2);
X      = S(3);
X_star = S(4);


dZ1dt     = mu - (eta * Z1 * Z2) - (delta * Z1);
dZ2dt     = (theta * X) - (eta * Z1 * Z2) - (delta * Z2);
dXdt      = (k * Z1) - ((gamma + delta) * X) - (alpha_1 * X) + (alpha_2 * X_star);
dX_stardt = (alpha_1 * X) - (alpha_2 * X_star) - (gamma_star * X_star);


dS = [dZ1dt; dZ2dt; dXdt; dX_stardt];

end 