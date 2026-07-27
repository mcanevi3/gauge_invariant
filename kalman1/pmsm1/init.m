% init.m
% Initialization script for Armature-Controlled DC Motor
% System: dx/dt = A*x + Bu*u + Bw*w, y = C*x + Du*u + Dw*w
% States (x): [Armature Current (i_a); Rotor Speed (w_m)]
% Control (u): Armature Voltage (v_a)
% Disturbance (w): Load Torque (T_L)
% Output (y): Measured Armature Current (i_a)

clear; clc;

%% 1. Physical Parameters (Based on 24V Brushed DC Servo Motor)
% Electrical Parameters
Ra = 1.2;       % Armature resistance [Ohms]
La = 0.0015;    % Armature inductance [Henries] (1.5 mH)
Kt = 0.05;      % Motor torque constant [Nm/A]
Ke = 0.05;      % Back-EMF constant [V/(rad/s)] (Equals Kt in SI units)

% Mechanical Parameters
J = 1.5e-4;     % Rotor moment of inertia [kg*m^2]
B = 1.0e-5;     % Viscous friction coefficient [Nms/rad]

%% 2. State-Space Matrices
% Continuous-time LTI system matrices
A = [-Ra/La,  -Ke/La;
      Kt/J,    -B/J ];

Bu = [ 1/La;
       0   ];

Bw = [ 0;
      -1/J ];

C  = [ 1, 0 ];  % Only current is measured

Du = 0;
Dw = 0;

%% 3. Create MATLAB State-Space Objects
% Model mapping control input (Voltage) to measured output (Current)
sys_u = ss(A, Bu, C, Du);
sys_u.StateName  = {'i_a', 'w_m'};
sys_u.InputName  = {'v_a'};
sys_u.OutputName = {'i_a_meas'};

% Model mapping disturbance (Load Torque) to measured output (Current)
sys_w = ss(A, Bw, C, Dw);
sys_w.StateName  = {'i_a', 'w_m'};
sys_w.InputName  = {'T_L'};
sys_w.OutputName = {'i_a_meas'};

%% 4. Display System Information
disp('--------------------------------------------------');
disp('DC Motor Initialization Complete.');
disp('System matrices A, Bu, Bw, C, Du, Dw loaded.');
disp(' ');
fprintf('Electrical Time Constant (tau_e): %.3f ms\n', (La/Ra)*1000);
fprintf('Mechanical Time Constant (tau_m): %.1f ms\n', (J*Ra/(Kt*Ke))*1000);
disp('--------------------------------------------------');