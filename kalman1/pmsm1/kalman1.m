% Augmented_KF_DCMotor.m
% Discrete Kalman Filter for State and Disturbance Estimation
clear; clc; close all;

%% 1. Motor Parameters (from init.m)
Ra = 1.2;       % Armature resistance [Ohms]
La = 0.0015;    % Armature inductance [H]
Kt = 0.05;      % Torque constant [Nm/A]
Ke = 0.05;      % Back-EMF constant [V/(rad/s)]
J  = 1.5e-4;    % Rotor inertia [kg*m^2]
B  = 1.0e-5;    % Viscous friction [Nms/rad]

%% 2. Augmented Continuous-Time System
% x_aug = [i_a; w_m; T_L]
A_aug = [-Ra/La, -Ke/La,    0;
          Kt/J,   -B/J,  -1/J;
          0,       0,       0];

Bu_aug = [1/La; 
          0; 
          0];

C_aug = [1, 0, 0];  % Only i_a is measured
D_aug = 0;

%% 3. Discretization
Ts = 1e-3; % 1 ms sampling time
sys_c = ss(A_aug, Bu_aug, C_aug, D_aug);
sys_d = c2d(sys_c, Ts, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

%% 4. Kalman Filter Initialization
% Tuning Matrices
% Q represents confidence in our model. High value on T_L (index 3,3) 
% tells the filter to expect the disturbance to change.
Q = diag([1e-4, 1e-4, 5e-2]); 

% R represents measurement noise covariance (sensor noise variance)
R = 1e-3; 

% Initial Conditions
x_est = [0; 0; 0];    % Initial state estimate
P     = eye(3);       % Initial error covariance

%% 5. Simulation Setup
t_sim = 4;                  % Simulation time [s]
N     = t_sim / Ts;         % Number of samples
t     = (0:N-1) * Ts;       % Time vector

% Inputs and Disturbances
v_a = 24 * ones(1, N);      % 24V step input
T_L_true = zeros(1, N);
T_L_true(t >= 0.5) = 1;   % 0.2 Nm load step at t = 0.5s
T_L_true(t >= 1.5) = -1;   % 0.2 Nm load step at t = 0.5s
T_L_true(t >= 2.5) = 0;   % 0.2 Nm load step at t = 0.5s

% True plant state and measurement initialization
x_true = [0; 0; 0];
y_meas = zeros(1, N);

% Storage arrays for plotting
x_true_history = zeros(3, N);
x_est_history  = zeros(3, N);

%% 6. Simulation and Kalman Filter Loop
for k = 1:N
    % --- PLANT SIMULATION ---
    % Add simulated process noise and measurement noise
    w_proc = [1e-3*randn; 1e-3*randn; 0]; % Noise on true states
    v_meas = sqrt(R) * randn;             % Sensor noise
    
    % Force the true disturbance into the state vector for simulation
    x_true(3) = T_L_true(k); 
    
    % Output measurement
    y_meas(k) = Cd * x_true + v_meas;
    
    % Store true states
    x_true_history(:, k) = x_true;
    
    % Plant state evolution (using discrete equations)
    x_true = Ad * x_true + Bd * v_a(k) + w_proc;
    
    % --- KALMAN FILTER ---
    % 1. Predict Step (Time Update)
    x_pred = Ad * x_est + Bd * v_a(k);
    P_pred = Ad * P * Ad' + Q;
    
    % 2. Update Step (Measurement Update)
    K = P_pred * Cd' / (Cd * P_pred * Cd' + R);    % Kalman Gain
    x_est = x_pred + K * (y_meas(k) - Cd * x_pred); % State Update
    P = (eye(3) - K * Cd) * P_pred;                 % Covariance Update
    
    % Store estimates
    x_est_history(:, k) = x_est;
end

%% 7. Plotting Results
figure(1);

subplot(3,1,1);
plot(t, x_true_history(1,:), 'b', t, x_est_history(1,:), 'r--', 'LineWidth', 1.5);
ylabel('Current (A)'); title('Armature Current (Measured vs Estimated)');
legend('True', 'Estimated'); grid on;

subplot(3,1,2);
plot(t, x_true_history(2,:), 'b', t, x_est_history(2,:), 'r--', 'LineWidth', 1.5);
ylabel('Speed (rad/s)'); title('Rotor Speed (Estimated)');
legend('True', 'Estimated'); grid on;

subplot(3,1,3);
plot(t, T_L_true, 'b', t, x_est_history(3,:), 'r--', 'LineWidth', 1.5);
ylabel('Torque (Nm)'); xlabel('Time (s)'); title('Load Torque Disturbance (Estimated)');
legend('True', 'Estimated'); grid on;