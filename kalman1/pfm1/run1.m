% run_filter_w_only.m
clear; clc; close all;

% 1. Load System Matrices
init; 

% 2. Simulation Setup
N = 1000;              % Number of simulation steps
t = (0:N-1) * tau;     % Time vector
u = sin(2*pi*5*t);     % Example control input

% 3. Noise Covariances
Q_x = diag([1e-6, 1e-6]); 
Q_w = 1e-4;               % Process noise covariance for disturbance
R   = diag([1e-4, 1e-4, 1e-4]); 

% 4. Projection Matrix Setup
Pi = eye(n_y) - C * pinv(C'*C) * C';
Pi_Dw = Pi * Dw;
R_bar = Pi * R * Pi';     % Projected measurement noise covariance

% 5. Full Pre-allocation and Signal Generation
x_true = zeros(n_x, N);
y_meas = zeros(n_y, N);
w_hat  = zeros(n_w, N);
w_true = zeros(n_w, N);

% Generate a complex disturbance profile offline
w_true(1, 201:400) = 2.0;                                     % Step
w_true(1, 401:700) = 2.0 + 1.5 * sin(2*pi*2*t(401:700));      % Sine wave offset
w_true(1, 701:N)   = linspace(2.0, -1.0, N-700);              % Ramp down

% Pre-generate noise arrays to avoid randn() overhead inside the loop
V_x = sqrt(diag(Q_x)) .* randn(n_x, N);
V_k = sqrt(diag(R)) .* randn(n_y, N);

% Initialize states and covariances
x_true(:,1) = [0; 0];
w_hat(:,1)  = 0;
P_w = eye(n_w);

% 6. Main Simulation and Filter Loop
for k = 2:N
    % --- SIMULATE TRUE SYSTEM ---
    % Notice how clean this is without signal generation inside the loop
    x_true(:,k) = Ad * x_true(:,k-1) + Bud * u(k-1) + Bwd * w_true(:,k-1) + V_x(:,k);
    y_meas(:,k) = C * x_true(:,k) + Dw * w_true(:,k) + V_k(:,k);
    
    
    % --- STAGE 1: DISTURBANCE ESTIMATION ONLY ---
    % Prediction
    w_pred = w_hat(:,k-1);
    P_w_pred = P_w + Q_w;
    
    % Projected measurement
    y_bar = Pi * y_meas(:,k);
    
    % Update
    S_w = Pi_Dw * P_w_pred * Pi_Dw' + R_bar;
    K_w = P_w_pred * Pi_Dw' * pinv(S_w); 
    
    w_hat(:,k) = w_pred + K_w * (y_bar - Pi_Dw * w_pred);
    P_w = (eye(n_w) - K_w * Pi_Dw) * P_w_pred;
end

% 7. Plot Results
figure('Name', 'Disturbance Estimation (Decoupled)');
plot(t, w_true(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, w_hat(1,:), 'r--', 'LineWidth', 1.5);
ylabel('w'); xlabel('Time (s)'); 
title('Disturbance: True vs Estimate (\Pi projection)');
legend('True Disturbance', 'Estimated Disturbance', 'Location', 'best');
grid on;